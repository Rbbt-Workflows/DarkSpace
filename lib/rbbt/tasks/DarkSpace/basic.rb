require 'rbbt/util/R'
module DarkSpace

  task :pmid_pair_by_resource => :tsv do

    dumper = TSV::Dumper.new :type => :double, :key_field => "PMID", :fields => ["IMEX", "BioGRID", "GO", "Reactome", "OmniPath_interactions", "OmniPath_ptm", "IID_pred", "STRING_pi", "STRING_tm", "EPMC", "EVEX"]
    dumper.init
    parser = TSV::Parser.new DarkSpace.triplets, :type => :double
    fields = parser.fields
    TSV.traverse parser, :into => dumper, :bar => self.progress_bar("Processing triplets") do |pair, values|
      pair = pair.first if Array === pair
      p1, p2 = pair.split("_")
      next if p1.nil? or p2.nil?

      pair = [p1, p2] * ":"

      res = [].extend MultipleResult
      Misc.zip_fields(values).each do |list|
        pmid, imex, reactome, tm_epmc, tm_pr_times_found, tm_dm_times_found, evex, evex_score, biogrid, go_ipi, omnipath_interactions, omnipath_ptm, iid_pred, string_pi_score_ave, string_pi, string_tm_score_ave, string_textmining, noncur_prot, pmid2 = list
        next unless pmid =~ /^\d{6,}$/ 
        triplet_values = [imex, biogrid, go_ipi, reactome, omnipath_interactions, omnipath_ptm, iid_pred, string_pi, string_textmining, tm_epmc, evex].collect{|_v| _v == "1" ? [pair] : []}
        res << [pmid, triplet_values]
      end
      res
    end

    io = TSV.collapse_stream(dumper.stream).stream
    Misc.open_pipe do |sin|
      while line = io.gets
        parts = line.split("\t").collect{|p| p.split("|").compact.uniq.reject{|s| s.empty?} * "|" }
        sin.puts parts * "\t"
      end
    end
  end

  task :pmid_tm_scores => :tsv do
    dumper = TSV::Dumper.new :key_field => "PMID", :fields => ["EPMC PR score", "EPMC DR score", "EVEX score", "STRING score"], :type => :double, :namespace => DarkSpace.organism, :cast => :to_f
    dumper.init
    parser = TSV::Parser.new DarkSpace.triplets, :type => :double, :fields => ["pmid", "tm_pr_times_found", "tm_dm_times_found", "evex_score", "STRING_tm_score_ave"]
    fields = parser.fields
    TSV.traverse parser, :into => dumper, :bar => self.progress_bar("Processing triplets") do |pair, values|
      pair = pair.first if Array === pair
      p1, p2 = pair.split("_")
      next if p1.nil? or p2.nil?

      res = []
      Misc.zip_fields(values).each do |pmid, tm_pr, tm_dm, evex, string|
        next unless pmid =~ /^\d{6,}$/
        res << [pmid, [tm_pr, tm_dm, evex, string].collect{|v| v == "0" ? nil : v.to_f }]
      end
      next if res.empty?
      res.extend MultipleResult
      res
    end

    TSV.collapse_stream(dumper.stream)
  end

  dep :pmid_pair_by_resource
  dep :pmid_tm_scores
  input :max, :boolean, "Use maximum of scores instead of average", true
  task :pmid_tm_relevance => :tsv do |max|

    require 'rbbt/util/R' 


    tsv = step(:pmid_pair_by_resource).load

    tsv.add_field "DB score" do |k,values|
      imex = values["IMEX"].any? ? 1 : 0
      ppi  = values.values_at("BioGRID", "GO").select{|l| l.any?}.length

      any_ppi = imex + ppi > 0 ? 1 : 0

      rest = values.values_at("Reactome", "OmniPath_interactions", "OmniPath_ptm", "IID_pred", "STRING_pi").select{|l| l.any?}.length

      imex * 10 + ppi * 10 + rest * 1
    end

    if max
      scores = step(:pmid_tm_scores).load.to_list{|values| values.empty? ? nil : Misc.max(values.collect{|v| v.to_f}) }
    else
      scores = step(:pmid_tm_scores).load.to_list{|values| values.empty? ? nil : Misc.mean(values.collect{|v| v.to_f}) }
    end

    tsv = tsv.select(scores.keys).attach(scores)

    pred = tsv.R <<-EOF
names(data) <- make.names(names(data))
data[is.na(data)] = 0
rbbt.require('randomForest')
d = data[sample(rownames(data), 10000),]
m = randomForest(DB.score ~ EPMC.PR.score + EPMC.DR.score + EVEX.score + STRING.score, d)
p = predict(m,data)
data = data.frame(p)
rownames(data) <- names(p)
names(data) <- c('Relevance')
    EOF

    pred.key_field = "PMID"
    tsv.attach pred
  end

  dep :predicted_ppi
  dep :pmid_pair_by_resource
  task :pmid_dark_space => :tsv do
    tsv = step(:pmid_pair_by_resource).load
    predicted_ppi = step(:predicted_ppi).load.collect{|p| p.split(":").sort.join(":")}

    tsv.monitor = true

    tsv.add_field "Known pairs" do |k, values|
      values.values_at("IMEX", "BioGRID", "GO").compact.flatten
    end

    tsv.add_field "Pathway pairs" do |k, values|
      values.values_at("Reactome", "OmniPath_interactions", "OmniPath_ptm", "IID_pred", "STRING_pi").compact.flatten
    end

    tsv.add_field "Text-Mining pairs" do |k, values|
      values.values_at("EPMC", "EVEX", "STRING_textmining").compact.flatten
    end

    tsv.add_field "Dark Space" do |k,values|
      values["Text-Mining pairs"] - values["Known pairs"]
    end

    uni2name = Organism.identifiers(DarkSpace.organism).index(:target => "Associated Gene Name", :fields => "UniProt/SwissProt Accession", :order => true, :persist => true)
    tsv.add_field "Predicted Dark Space" do |k,values|
      values["Dark Space"].select do |pair|
        parts = pair.split(":").collect{|p| uni2name[p]}
        next false if parts.compact.length < 2
        name_pair = parts.sort.join(":")
        predicted_ppi.include? name_pair
      end
    end

    tsv.add_field "Pathway matched Dark Space" do |k,values|
      proteins = [values["Known pairs"], values["Pathway pairs"]].compact.flatten.collect{|pair| pair.split(":") }.flatten.uniq
      values["Dark Space"].select do |pair|
        p1, p2 = pair.split(":")
        proteins.include?(p1) && proteins.include?(p2)
      end
    end

    tsv.add_field "Pathway partial-matched Dark Space" do |k,values|
      proteins = [values["Known pairs"], values["Pathway pairs"]].compact.flatten.collect{|pair| pair.split(":") }.flatten.uniq
      values["Dark Space"].select do |pair|
        p1, p2 = pair.split(":")
        proteins.include?(p1) || proteins.include?(p2)
      end
    end

    tsv
  end

  dep :pmid_pair_by_resource
  task :protein_counts => :tsv do
    tsv = step(:pmid_pair_by_resource).load 

    protein_counts = TSV.setup({}, :key_field => "UniProt/SwissProt Accession", :fields => ["Counts"], :type => :single, :cast => :to_i)
    tsv.through do |pmid,values|
      values.values_at("IMEX", "BioGRID", "GO").flatten.compact.uniq.collect{|pair| pair.split(":")}.flatten.each do |protein|
        protein_counts[protein] ||= 0
        protein_counts[protein] += 1
      end
    end

    protein_counts
  end

  input :score, :float, "Minumum PPI prediction score", 0
  task :predicted_ppi => :array do |score|
    tsv = PIPS.data.tsv
    tsv.unzip.select("Score"){|s| s.to_f >= score}.keys
  end


  dep :pmid_dark_space, :compute => :produce
  dep :pmid_tm_relevance, :compute => :produce
  dep :protein_counts, :compute => :produce
  task :pmid_ranks => :tsv do
    tsv = step(:pmid_dark_space).load

    tsv.monitor = true
    tsv.attach step(:pmid_tm_relevance)

    protein_counts = step(:protein_counts).load

    tsv.add_field "Dark Space interest" do |k,values|
      proteins = values["Pathway matched Dark Space"].collect{|pair| pair.split(":")}
      proteins.inject(0){|acc,protein| acc += 1.0 / ((protein_counts[protein] || 0) + 1) }
    end

    tsv.add_field "Dark Space (partial) interest" do |k,values|
      proteins = values["Pathway partial-matched Dark Space"].collect{|pair| pair.split(":")}
      proteins.inject(0){|acc,protein| acc += 1.0 / ((protein_counts[protein] || 0) + 1) }
    end

    tsv.add_field "Predicted Dark Space interest" do |k,values|
      proteins = values["Predicted Dark Space"].collect{|pair| pair.split(":")}
      proteins.inject(0){|acc,protein| acc += 1.0 / ((protein_counts[protein] || 0) + 1) }
    end

    tsv
  end

  task :pmid_validation => :tsv do
    tsv1 = Rbbt.manual_evaluation.round1["eval_miguel_rd1_validation.txt"].tsv :type => :list, :cast => :to_s, :fields => ["contains interactions?"]
    tsv1 = tsv1.add_field "Valid" do |k,v|
      v["contains interactions?"] == "yes"
    end.slice("Valid")

    tsv2 = Rbbt.manual_evaluation.round2["eval_miguel_rd2_validation.txt"].tsv :type => :list, :cast => :to_s, :fields => ["contains interactions?"]
    tsv2 = tsv2.add_field "Valid" do |k,v|
      v["contains interactions?"] == "yes"
    end.slice("Valid")

    tsv3 = Rbbt.manual_evaluation.round3["eval_miguel_rd3_validation.txt"].tsv :type => :list, :cast => :to_s, :fields => ["contains_interactions"]
    tsv3 = tsv3.add_field "Valid" do |k,v|
      v["contains_interactions"] == "yes"
    end.slice("Valid")

    tsv1.add_field "Source" do 
      "Round1"
    end
    tsv2.add_field "Source" do 
      "Round2"
    end
    tsv3.add_field "Source" do 
      "Round3"
    end

    tsv1.merge! tsv2
    tsv1.merge! tsv3

    tsv1
  end

  dep :pmid_validation
  dep :pmid_ranks
  task :eval => :tsv do
    ranks = step(:pmid_ranks).load
    validation = step(:pmid_validation).load

    valid = validation.select("Valid" => "true").keys
    nonvalid = validation.select("Valid" => "false").keys

    valid_relevance = ranks.values_at(*valid).collect{|v| v["Relevance"].first.to_f}
    nonvalid_relevance = ranks.values_at(*nonvalid).collect{|v| v["Relevance"].first.to_f}

    require 'rbbt/util/R'
    tsv = TSV.setup({}, :key_field => "Statistic", :fields => ["Value"], :type => :single)
    tsv["Valid mean"] = Misc.mean valid_relevance
    tsv["NonValid mean"] = Misc.mean nonvalid_relevance
    tsv["p-value"] = R.eval("t.test(#{R.ruby2R valid_relevance}, #{R.ruby2R nonvalid_relevance})$p.value")
    tsv
  end

  dep :pmid_validation
  dep :pmid_ranks
  extension :png
  task :eval_ROC => :binary do
    ranks = step(:pmid_ranks).load
    validation = step(:pmid_validation).load

    valid = validation.select("Valid" => "true").keys
    nonvalid = validation.select("Valid" => "false").keys

    valid_relevance = ranks.values_at(*valid).collect{|v| v["Relevance"].first.to_f}
    nonvalid_relevance = ranks.values_at(*nonvalid).collect{|v| v["Relevance"].first.to_f}

    all_scores = (valid_relevance + nonvalid_relevance).flatten.sort
    score_values = {}
    all_scores.each do |score|
      tp = valid_relevance.select{|s| s >= score}.length
      fp = nonvalid_relevance.select{|s| s >= score}.length

      tn = nonvalid_relevance.select{|s| s < score}.length
      fn = valid_relevance.select{|s| s < score}.length

      precision = tp.to_f / (tp + fp)
      recall = tpr = tp.to_f / valid_relevance.length

      fpr = fp.to_f / nonvalid_relevance.length

      score_values[score] = [tp, fp, tn, fn, tpr, fpr, precision, recall]
    end

    recall = score_values.sort_by{|score,values| score}.collect{|score,values| values[-1] }.flatten
    precision = score_values.sort_by{|score,values| score}.collect{|score,values| values[-2] }.flatten
    fpr = score_values.sort_by{|score,values| score}.collect{|score,values| values[-3] }.flatten
    tpr = score_values.sort_by{|score,values| score}.collect{|score,values| values[-4] }.flatten

    R.run <<-EOF, [:svg]

precision = #{R.ruby2R precision}
recall = #{R.ruby2R recall}
tpr = #{R.ruby2R tpr}
fpr = #{R.ruby2R fpr}
png('#{self.path}', width=400, height=800)
par(mfrow=c(2,1))
plot(recall, precision, xlim=c(0,1), ylim=c(0,1))
plot(fpr, tpr, xlim=c(0,1), ylim=c(0,1))
dev.off
    EOF
    nil
  end

  dep :pmid_tm_relevance
  extension :svg
  task :pmid_tm_relevance_plot => :tsv do

    tsv = step(:pmid_tm_relevance).load

    tsv.R <<-EOF, [:svg]
rbbt.require('ggplot2')
names(data) <- make.names(names(data))
p = ggplot(data) + geom_smooth(aes(x=Relevance, y=DB.score)) + geom_point(aes(x=Relevance, y=DB.score), alpha=0.1)

rbbt.SVG.save.fast('#{self.path}', p, width=5, height=5)
data = NULL
    EOF

    nil
  end


end

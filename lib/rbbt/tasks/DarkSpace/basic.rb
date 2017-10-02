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
  task :pmid_tm_relevance => :tsv do

    require 'rbbt/util/R' 


    tsv = step(:pmid_pair_by_resource).load

    tsv.add_field "DB score" do |k,values|
      ppi  = values.values_at("IMEX", "BioGRID", "GO").select{|l| l.any?}.length
      rest = values.values_at("Reactome", "OmniPath_interactions", "OmniPath_ptm", "IID_pred", "STRING_pi").select{|l| l.any?}.length
      ppi * 10 + rest
    end

    scores = step(:pmid_tm_scores).load.to_list{|values| values.empty? ? nil : Misc.mean(values.collect{|v| v.to_f}) }

    tsv = tsv.attach(scores)

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

  dep :pmid_pair_by_resource
  task :pmid_dark_space => :tsv do
    tsv = step(:pmid_pair_by_resource).load

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

    tsv
  end


end

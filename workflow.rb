
require 'rbbt-util'
require 'rbbt/workflow'

Misc.add_libdir if __FILE__ == $0

require 'rbbt/sources/DarkSpace'

module DarkSpace
  extend Workflow

  task :protein_evidence => :tsv do
    dumper = TSV::Dumper.new :key_field => "UniProt/SwissProt Accession", :fields => ["Partner", "Resource"], :type => :double, :namespace => DarkSpace.organism
    dumper.init
    parser = TSV::Parser.new DarkSpace.pairs, :type => :list
    fields = parser.fields
    TSV.traverse parser, :into => dumper, :bar => self.progress_bar("Processing pairs") do |pair, values|
      pair = pair.first if Array === pair
      p1, p2 = pair.split("_")
      next if p1.nil? or p2.nil?
      res = []
      fields.zip(values).each do |field, value|
        next if field.downcase.include?('evex') or field.include?('textmining') or field.include?('tm_')
        next if value == '0'
        res << [p1, [p2, field]]
        res << [p2, [p1, field]] unless p1 == p2
      end
      next if res.empty?
      res.extend MultipleResult
      res
    end
    TSV.collapse_stream(dumper.stream)
  end

  task :pmid_proteins => :tsv do
    dumper = TSV::Dumper.new :key_field => "PMID", :fields => ["UniProt/SwissProt Accession"], :type => :flat, :namespace => DarkSpace.organism
    dumper.init
    parser = TSV::Parser.new DarkSpace.triplets, :type => :list
    fields = parser.fields
    TSV.traverse parser, :into => dumper, :bar => self.progress_bar("Processing triplets") do |pair, values|
      pair = pair.first if Array === pair
      p1, p2 = pair.split("_")
      next if p1.nil? or p2.nil?
      pmids = values.first.split("|")
      res = []

      pmids.each do |pmid|
        res << [pmid, [p1, p2].uniq * "|"]
      end

      res.extend MultipleResult
      res
    end
    TSV.collapse_stream(dumper.stream)
  end

  task :pmid_coverage => :tsv do
    dumper = TSV::Dumper.new :key_field => "UniProt/SwissProt Accession", :fields => ["Resource"], :type => :flat, :namespace => DarkSpace.organism
    dumper.init
    parser = TSV::Parser.new DarkSpace.pmids, :type => :list
    fields = parser.fields
    TSV.traverse parser, :into => dumper, :bar => self.progress_bar("Processing pairs") do |pmid, values|
      pmid = pmid.first if Array === pmid
      res = []
      fields.zip(values).each do |field, value|
        next if field.downcase.include?('evex') or field.include?('textmining') or field.include?('tm_')
        next if value == '0'
        res << [pmid, [field]]
      end
      next if res.empty?
      res.extend MultipleResult
      res
    end
    TSV.collapse_stream(dumper.stream)
  end

  task :pmid_tm_scores => :tsv do
    dumper = TSV::Dumper.new :key_field => "PMID", :fields => ["EPMC PR", "EPMC DR", "EVEX", "STRING"], :type => :double, :namespace => DarkSpace.organism, :cast => :to_f
    dumper.init
    parser = TSV::Parser.new DarkSpace.triplets, :type => :double, :fields => ["pmid", "tm_pr_times_found", "tm_dm_times_found", "evex_score", "STRING_tm_score_ave"]
    fields = parser.fields
    TSV.traverse parser, :into => dumper, :bar => self.progress_bar("Processing triplets") do |pair, values|
      pair = pair.first if Array === pair
      p1, p2 = pair.split("_")
      next if p1.nil? or p2.nil?

      res = []
      Misc.zip_fields(values).each do |pmid, tm_pr, tm_dm, evex, string|
        res << [pmid, [tm_pr, tm_dm, evex, string].collect{|v| v == "0" ? nil : v.to_f }]
      end
      next if res.empty?
      res.extend MultipleResult
      res
    end

    TSV.collapse_stream(dumper.stream)
  end

  dep :pmid_coverage, :compute => :produce
  dep :pmid_tm_scores
  task :pmid_tm_scores_db => :tsv do
    pmid_coverage = step(:pmid_coverage).load

    pmid_tm_scores = step(:pmid_tm_scores).load

    pmid_tm_scores = pmid_tm_scores.to_list{|v| Misc.mean(v.collect{|v| v.nil? ? nil : v.to_f })}


    pmid_tm_scores.monitor = true
    pmid_tm_scores.add_field "IMEX" do |pmid,values|
      pmid_coverage.include?(pmid) && pmid_coverage[pmid].include?('imex')
    end

    pmid_tm_scores.add_field "Num. DB" do |pmid,values|
      pmid_coverage.include?(pmid) ? pmid_coverage[pmid].length : 0
    end

    pmid_tm_scores
  end

  dep :pmid_tm_scores_db
  task :pmid_tm_relevance => :tsv do

    require 'rbbt/util/R' 

    tsv = step(:pmid_tm_scores_db).load
    tsv.add_field "DB Score" do |k,values|
      values["IMEX"] == 'true' ? (10 + values["Num. DB"].to_f) : values["Num. DB"].to_f
    end

    pred = tsv.R <<-EOF
names(data) <- make.names(names(data))
data[is.na(data)] = 0
rbbt.require('randomForest')
d = data[sample(rownames(data), 10000),]
m = randomForest(DB.Score ~ EPMC.PR + EPMC.DR + EVEX + STRING, d)
p = predict(m,data)
data = data.frame(p)
rownames(data) <- names(p)
names(data) <- c('Relevance')
    EOF

    pred
  end

  dep :pmid_coverage, :compute => :produce
  dep :pmid_tm_relevance
  dep :pmid_proteins
  dep :protein_evidence
  task :pmid_rank => :tsv do
    pmid_coverage = step(:pmid_coverage).load
    scores = step(:pmid_tm_relevance).load
    pmid_proteins = step(:pmid_proteins).load
    protein_evidence = step(:protein_evidence).load

    scores.monitor = true

    scores.add_field "IMEX" do |pmid,values|
      pmid_coverage.include?(pmid) ? pmid_coverage[pmid].include?("imex") : false
    end

    scores.add_field "Coverage" do |pmid,values|
      pmid_coverage.include?(pmid) ? pmid_coverage[pmid].length : 0
    end

    scores.add_field "Proteins" do |pmid,values|
      pmid_proteins[pmid].uniq
    end

    scores.add_field "Protein interest" do |pmid,values|
      interest = values["Proteins"].collect{|protein|
        evidence = protein_evidence.include?(protein) ? protein_evidence[protein].last : []

        1 / (1 + evidence.length)
      }
      Misc.sum(interest)
    end

    scores.add_field "Combined score" do |pmid,values|
      if values["IMEX"].to_s == "true"
        0
      else
        i = values["Protein interest"].to_f
        p = values["Relevance"].to_f
        c = values["Coverage"].to_f + 1
        i * p / c
      end
    end

    scores = scores.select(scores.keys.select{|k| k =~ /^\d{6,10}$/})

    scores
  end

end

#require 'DarkSpace/tasks/basic.rb'

#require 'rbbt/knowledge_base/DarkSpace'
#require 'rbbt/entity/DarkSpace'


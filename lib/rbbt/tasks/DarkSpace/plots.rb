module DarkSpace

#  dep :pathway_pmid_mentions
#  dep :pmid_coverage
#  extension :svg
#  task :pathway_pmid_db_coverage => :text do
#    pmids = step(:pathway_pmid_mentions).load.keys
#    all_pmids = DarkSpace.pmids.tsv.keys.select{|p| p =~ /^\d{6,}$/}
#    coverage = step(:pmid_coverage).load
#
#    coverage.delete "-"
#
#    coverage_melted = coverage.annotate({})
#    coverage_melted.fields = coverage.all_fields
#
#    found = Set.new
#    coverage.each do |pmid,dbs|
#      found << pmid
#      dbs.each do |db|
#        value = [pmid, db]
#        id = Misc.obj2digest value
#        coverage_melted[id] = value
#      end
#    end
#
#    all_pmids.each do |pmid|
#      next if found.include? pmid
#      value = [pmid, "NONE"]
#      id = Misc.obj2digest value
#      coverage_melted[id] = value
#    end
#
#
#    require 'rbbt/util/R'
#
#    TmpFile.with_file(pmids * "\n") do |pmid_file|
#      coverage_melted.R_interactive <<-EOF
#pmids = scan('#{pmid_file}')
#coverage.m = data
#      EOF
#    end
#  end
end

require 'rbbt-util'
require 'rbbt/resource'
require 'rbbt/sources/organism'

module DarkSpace
  extend Resource
  self.subdir = 'share/databases/DarkSpace'

  def self.organism(org="Hsa")
    Organism.default_code(org)
  end

  #self.search_paths = {}
  #self.search_paths[:default] = :lib


  DarkSpace.claim DarkSpace.pairs, :proc do 
    Rbbt.modules.darkspaceproject.dsp_comparison.results["paircomp_table_final.txt"].tsv :header_hash => ""
  end

  DarkSpace.claim DarkSpace.pmids, :proc do 
    Rbbt.modules.darkspaceproject.dsp_comparison.results["pubcomp_table_final.txt"].tsv :header_hash => ""
  end

  DarkSpace.claim DarkSpace.triplets, :proc do 
    tsv = TmpFile.with_file do |dir|
      Path.setup(dir)
      tar = Rbbt.modules.darkspaceproject.dsp_comparison.results["pubpaircomp_table_final.txt.tar.gz"]
      Misc.untar(tar.open, dir)
      dir.glob("*.txt").first.tsv :header_hash => "", :merge => true
    end
    tsv.zip("pmid", true)
  end
end

iif DarkSpace.pairs.produce.find if __FILE__ == $0
iif DarkSpace.pmids.produce.find if __FILE__ == $0
iif DarkSpace.triplets.produce.find if __FILE__ == $0


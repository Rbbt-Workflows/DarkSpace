
require 'rbbt-util'
require 'rbbt/resource'

module PIPS
  extend Resource
  self.subdir = 'share/databases/PIPS'

  #def self.organism(org="Hsa")
  #  Organism.default_code(org)
  #end

  #self.search_paths = {}
  #self.search_paths[:default] = :lib


  PIPS.claim PIPS.data, :proc do 
    url = "http://www.compbio.dundee.ac.uk/www-pips/downloads/PredictedInteractions100.txt"
    tsv = TSV.open(url, :header_hash => "", :merge => true)
    tsv.key_field = "Protein 1 (Associated Gene Name)"
    tsv.fields = ["ACC1", "Protein 2 (Associated Gene Name)", "ACC2", "Score"]
    tsv = tsv.slice ["Protein 2 (Associated Gene Name)", "Score"]
    tsv
  end
end

iif PIPS.data.produce(true).find if __FILE__ == $0


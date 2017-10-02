
require 'rbbt-util'
require 'rbbt/workflow'

Misc.add_libdir if __FILE__ == $0

require 'rbbt/sources/DarkSpace'

module DarkSpace
  extend Workflow
end

require 'rbbt/tasks/DarkSpace/basic.rb'
require 'rbbt/tasks/DarkSpace/plots.rb'

#require 'rbbt/knowledge_base/DarkSpace'
#require 'rbbt/entity/DarkSpace'



Workflow.require_workflow "Genomics"

get '/pmid/info/:pmid' do
  a = PubMed.get_article(params[:pmid])
  info = {:title => a.title, :abstract => a.abstract, :year => a.year}

  content_type :json
  halt 200, info.to_json
end

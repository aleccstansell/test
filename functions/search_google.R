install.packages("searcher")
library(searcher)
gene_list

for (1:length(gene_list)) {
  search_site(query, site = c("google"), rlang = TRUE)
}



?search_site


search_site(query, site = c("google", "bing", "stackoverflow", "so", "github",
                            "gh", "duckduckgo", "ddg", "bitbucket", "bb", "ixquick"), rlang = TRUE)


?paste

term <- paste(gene_list[1,], "HIV", sep = " ", collapse)



for (i in 150:210) {
  term <- paste(gene_list[i,], "HIV-1", sep = " ", collapse = NULL)
  search_site(term, "google", rlang = FALSE)
  
}

search_site(term, "google", rlang = FALSE)








search_google(gene_list[1,])


gene_list[1,]

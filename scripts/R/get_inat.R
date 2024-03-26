# rejiggering source from rinat
# jby 2021.12.03


# https://api.inaturalist.org/v1/observations?quality_grade=research&place_id=53170&taxon_id=47785&term_id=12&term_value_id=13

# term id: 12 for "Plant Phenology" then term_id_value: 13 =Flowering, 14 =Fruiting, 15 =Flower Budding

library(rinat)

get_inat_obs <- function (query = NULL, taxon_name = NULL, taxon_id = NULL, place_id = NULL,  quality = NULL, geo = NULL, year = NULL, month = NULL, day = NULL, term_id = NULL, term_value_id = NULL, without_term_value_id = NULL, bounds = NULL, maxresults = 100, meta = FALSE) 
{
    if (!curl::has_internet()) {
        message("No Internet connection.")
        return(invisible(NULL))
    }
    base_url <- "http://www.inaturalist.org/"
    if (httr::http_error(base_url)) {
        message("iNaturalist API is unavailable.")
        return(invisible(NULL))
    }
    arg_list <- list(query, taxon_name, taxon_id, place_id, quality, 
        geo, year, month, day, bounds)
    arg_vals <- lapply(arg_list, is.null)
    if (all(unlist(arg_vals))) {
        stop("All search parameters NULL. Please provide at least one.")
    }
    search <- ""
    if (!is.null(query)) {
        search <- paste0(search, "&q=", gsub(" ", "+", query))
    }
    if (!is.null(quality)) {
        if (!sum(grepl(quality, c("casual", "research")))) {
            stop("Please enter a valid quality flag,'casual' or 'research'.")
        }
        search <- paste0(search, "&quality_grade=", quality)
    }
    if (!is.null(taxon_name)) {
        search <- paste0(search, "&taxon_name=", gsub(" ", "+", 
            taxon_name))
    }
    if (!is.null(taxon_id)) {
        search <- paste0(search, "&taxon_id=", gsub(" ", "+", 
            taxon_id))
    }
	if (!is.null(term_id)) {
        search <- paste0(search, "&term_id=", term_id) # my new thing
    }
	if (!is.null(term_value_id)) {
        search <- paste0(search, "&term_value_id=", term_value_id) # my new thing
    }
	if (!is.null(without_term_value_id)) {
        search <- paste0(search, "&without_term_value_id=", without_term_value_id) # my new thing
    }

    if (!is.null(place_id)) {
        search <- paste0(search, "&place_id=", gsub(" ", "+", 
            place_id))
    }
    if (!is.null(geo) && geo) {
        search <- paste0(search, "&has[]=geo")
    }
    if (!is.null(year)) {
        if (length(year) > 1) {
            stop("You can only filter results by one year, please enter only one value for year.")
        }
        search <- paste0(search, "&year=", year)
    }
    if (!is.null(month)) {
        month <- as.numeric(month)
        if (is.na(month)) {
            stop("Please enter a month as a number between 1 and 12, not as a word.")
        }
        if (length(month) > 1) {
            stop("You can only filter results by one month, please enter only one value for month.")
        }
        if (month < 1 || month > 12) {
            stop("Please enter a valid month between 1 and 12")
        }
        search <- paste0(search, "&month=", month)
    }
    if (!is.null(day)) {
        day <- as.numeric(day)
        if (is.na(day)) {
            stop("Please enter a day as a number between 1 and 31, not as a word.")
        }
        if (length(day) > 1) {
            stop("You can only filter results by one day, please enter only one value for day.")
        }
        if (day < 1 || day > 31) {
            stop("Please enter a valid day between 1 and 31")
        }
        search <- paste0(search, "&day=", day)
    }
    if (!is.null(bounds)) {
        if (inherits(bounds, "sf")) {
            bounds_prep <- sf::st_bbox(bounds)
            bounds <- c(swlat = bounds_prep[2], swlng = bounds_prep[1], 
                nelat = bounds_prep[4], nelng = bounds_prep[3])
        }
        if (inherits(bounds, "Spatial")) {
            bounds_prep <- sp::bbox(bounds)
            bounds <- c(swlat = bounds_prep[2, 1], swlng = bounds_prep[1, 
                1], nelat = bounds_prep[2, 2], nelng = bounds_prep[1, 
                2])
        }
        if (length(bounds) != 4) {
            stop("Bounding box specifications must have 4 coordinates.")
        }
        search <- paste0(search, "&swlat=", bounds[1], "&swlng=", 
            bounds[2], "&nelat=", bounds[3], "&nelng=", bounds[4])
    }
    
    cat("Running search", search, "\n")
    
    if (maxresults > 10000) {
        stop("Please provide a maxresults value <= 10000.")
    }
    q_path <- "observations.csv"
    ping_path <- "observations.json"
    ping_query <- paste0(search, "&per_page=1&page=1")
    ping <- httr::GET(base_url, path = ping_path, query = ping_query)
    total_res <- as.numeric(ping$headers$`x-total-entries`)
    if (total_res == 0) {
        stop("Your search returned zero results. Either your species of interest has no records or you entered an invalid search.")
    }
    else if (total_res >= 2e+05) {
        stop("Your search returned too many results, please consider breaking it up into smaller chunks by year or month.")
    }
    else if (!is.null(bounds) && total_res >= 1e+05) {
        stop("Your search returned too many results, please consider breaking it up into smaller chunks by year or month.")
    }
    page_query <- paste0(search, "&per_page=200&page=1")
    data <- httr::GET(base_url, path = q_path, query = page_query)
    data <- inat_handle(data)
    data_out <- if (is.na(data)) 
        NA
    else read.csv(textConnection(data), stringsAsFactors = FALSE)
    if (total_res < maxresults) 
        maxresults <- total_res
    if (maxresults > 200) {
        for (i in 2:ceiling(maxresults/200)) {
            page_query <- paste0(search, "&per_page=200&page=", 
                i)
            data <- httr::GET(base_url, path = q_path, query = page_query)
            data <- inat_handle(data)
            data_out <- rbind(data_out, read.csv(textConnection(data), 
                stringsAsFactors = FALSE))
        }
    }
    if (is.data.frame(data_out)) {
        if (maxresults < dim(data_out)[1]) {
            data_out <- data_out[1:maxresults, ]
        }
    }
    if (meta) {
        return(list(meta = list(found = total_res, returned = nrow(data_out)), 
            data = data_out))
    }
    else {
        return(data_out)
    }
}

inat_handle <- function(x){
  res <- httr::content(x, as = "text")
  if(!x$headers$`content-type` == 'text/csv; charset=utf-8' || x$status_code > 202 || nchar(res) == 0 ){
    if(!x$headers$`content-type` == 'text/csv; charset=utf-8'){
      warning("Content type incorrect, should be 'text/csv; charset=utf-8'")
      NA
    }
    if(x$status_code > 202){
      warning(sprintf("Error: HTTP Status %s", x$status_code))
      NA
    }
    if(nchar(res) == 0){
      warning("No data found")
      NA
    }
  } else { res }
}


assignInNamespace("get_inat_obs", get_inat_obs, ns="rinat")

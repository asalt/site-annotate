# Initialize SQLite connection
initialize_cache_db <- function(db_path = "cache.sqlite") {
  if (!fs::file_exists(db_path)) {
    db <- DBI::dbConnect(RSQLite::SQLite(), db_path)
    DBI::dbExecute(db, "CREATE TABLE cache (hashval TEXT PRIMARY KEY, object BLOB)")
    DBI::dbDisconnect(db)
  }
  return(db_path)
}

write_to_cache <- function(object, hashval, db_path = "cache.sqlite") {
  db <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  on.exit(DBI::dbDisconnect(db))

  # Serialize the object
  serialized_obj <- serialize(object, NULL)

  # Insert or replace into the database
  DBI::dbExecute(
    db,
    "INSERT OR REPLACE INTO cache (hashval, object) VALUES (?, ?)",
    params = list(hashval, serialized_obj)
  )

  log_msg(msg = paste0("Saved object with hash: ", hashval, " to cache"))
}

load_from_cache <- function(hashval, db_path = "cache.sqlite") {
  db <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  on.exit(DBI::dbDisconnect(db))

  # Query the object
  result <- DBI::dbGetQuery(db, "SELECT object FROM cache WHERE hashval = ?", params = list(hashval))
  
  if (nrow(result) == 0) {
    log_msg(msg = paste0("Cache miss for hash: ", hashval))
    return(NULL)
  }

  log_msg(msg = paste0("Cache hit for hash: ", hashval))
  return(unserialize(result$object[[1]]))
}


get_hash_val <- function(obj) {
  rlang::hash(obj)
}

# hashval <- get_hash_val()
# 
# # Path to the cache database
# cache_db <- initialize_cache_db("cache.sqlite")
# 
# if (do_load) {
#   cached_result <- load_from_cache(hashval, db_path = cache_db)
#   if (!is.null(cached_result)) {
#     logger(msg = paste0("Cache hit: ", hashval))
#     return(cached_result)
#   } else {
#     return(NULL)
#   }
# } else if (save) {
#   write_to_cache(object = final_result, hashval = hashval, db_path = cache_db)
# }


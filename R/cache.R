suppressPackageStartupMessages(library(RSQLite))


`%cached_by%` <- function(compute_func, hashval) {

  compute_expr <- substitute(compute_func)

  # Check the cache
  cached_obj <- load_object_from_cache(hashval)

  if (!is.null(cached_obj)) {
    message(as.character(deparse(compute_expr)))
    message("Cache hit for: ", hashval)
    return(cached_obj)
  }

  # Perform computation
  message("Cache miss for hash: ", hashval)

  # Optional: Inspect or log the expression
  # can add some verbosity check here

  message("evaluating: ", deparse(compute_expr))

  result <- eval(compute_func)

  # Save to cache
  save_object_to_cache(result, hashval, notes = as.character(deparse(compute_expr)))

  return(result)
}



initialize_cache_db <- function(db_path = "cache.sqlite", close = TRUE){

  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = db_path)

  # Create tables if they don't exist
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS package_cache (
      package TEXT PRIMARY KEY,
      serialized_env BLOB,
      version TEXT,
      timestamp TEXT DEFAULT CURRENT_TIMESTAMP
    )
  ")

  # Create tables if they don't exist
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS object_cache (
      id INTEGER PRIMARY KEY AUTOINCREMENT,
      serialized_obj BLOB,
      object_hash TEXT NOT NULL UNIQUE,
      notes TEXT,
      timestamp TEXT DEFAULT CURRENT_TIMESTAMP
    )
  ")

  if (close) {
    DBI::dbDisconnect(con)
    message("Connection closed.")
    return(NULL)
  }
  return(con)
}


load_object_from_cache <- function(object_hash, db_path = "cache.sqlite", con = NULL, close = is.null(con)) {

  if (is.null(con)) {
    con <- DBI::dbConnect(RSQLite::SQLite(), dbname = db_path)
  }

  # Query the object cache
  result <- DBI::dbGetQuery(con, "
    SELECT serialized_obj FROM object_cache
    WHERE object_hash = ?
    LIMIT 1",
    params = list(object_hash)
  )

  if (close) {
    DBI::dbDisconnect(con)
    message("Connection closed.")
  }

  if (nrow(result) > 0) {
    # Deserialize and load the object
    serialized_obj <- result$serialized_obj[[1]]
    loaded_obj <- unserialize(serialized_obj)
    message("Object successfully loaded from cache.")
    return(loaded_obj)
  } else {
    message("No cached object found.")
    return(NULL)
  }
}


save_object_to_cache <- function(object, object_hash = NULL, notes = NULL, db_path = "cache.sqlite", con = NULL, close = is.null(con)) {

  if (is.null(con)) {
    con <- DBI::dbConnect(RSQLite::SQLite(), dbname = db_path)
  }

  if (is.null(notes)) {
    notes <- ""
  }

  if (length(notes) > 1) {
    notes <- paste(notes, collapse = "\n")
  }


  # Serialize the object
  serialized_obj <- serialize(object, NULL)
  serialized_obj_blob <- I(list(serialized_obj))

  # Get the object hash if not already provided
  if (is.null(object_hash)) {
    object_hash <- digest::digest(serialized_obj)
  }

  # Insert or update the object cache
  DBI::dbExecute(
    con,
    "INSERT INTO object_cache (serialized_obj, object_hash, notes, timestamp)
     VALUES (?, ?, ?, CURRENT_TIMESTAMP)
     ON CONFLICT(object_hash) DO UPDATE SET
       serialized_obj = excluded.serialized_obj,
       timestamp = CURRENT_TIMESTAMP",
    params = list(serialized_obj_blob, object_hash, notes)
  )

  if (close) {
    DBI::dbDisconnect(con)
    message("Connection closed.")
  }
}

save_package_to_db <- function(package, db_path = "cache.sqlite", con = NULL, close = is.null(con)) {

  if (is.null(con)) {
    con <- DBI::dbConnect(RSQLite::SQLite(), dbname = db_path)
  }

  # Load the package
  suppressPackageStartupMessages(library(package, character.only = TRUE))

  # Serialize the environment of the package
  package_env <- as.environment(paste0("package:", package))
  serialized_env <- I(serialize(as.list(package_env), NULL))
  serialized_env_blob <- I(list(serialized_env))

  # Get the package version
  package_version <- as.character(packageVersion(package))

  # Insert or update the package cache
  DBI::dbExecute(
    con,
    "INSERT INTO package_cache (package, serialized_env, version, timestamp)
     VALUES (?, ?, ?, CURRENT_TIMESTAMP)
     ON CONFLICT(package) DO UPDATE SET
       serialized_env = excluded.serialized_env,
       version = excluded.version,
       timestamp = CURRENT_TIMESTAMP",
    params = list(package, serialized_env_blob, package_version)
  )

  if (close) {
    DBI::dbDisconnect(con)
  }

  # DBI::dbDisconnect(con)
}



log_package_versions_to_db <- function(db_path = "cache.sqlite", con = NULL) {

  if (is.null(con)) {
    con <- DBI::dbConnect(RSQLite::SQLite(), dbname = db_path)
  }
  # Get installed package versions
  installed_versions <- data.frame(
    package = rownames(installed.packages()),
    version = as.character(installed.packages()[, "Version"]),
    stringsAsFactors = FALSE
  )

  # Update or insert package versions into the database
  DBI::dbWriteTable(con, "package_versions", installed_versions, append = TRUE, overwrite = FALSE)

  DBI::dbDisconnect(con)
}

load_package_from_db <- function(package, db_path = "cache.sqlite", con = NULL, close = is.null(con)) {

  close <- is.null(con)
  if (is.null(con)) {
    con <- DBI::dbConnect(RSQLite::SQLite(), dbname = db_path)
  }

  # Query the package cache
  result <- DBI::dbGetQuery(con, "
    SELECT serialized_env FROM package_cache
    WHERE package = ?
    LIMIT 1",
    params = list(package)
  )

  if (close) {
    DBI::dbDisconnect(con)
    message("Connection closed.")
  }


  if (nrow(result) > 0) {

    # Deserialize and load the package environment
    serialized_env <- result$serialized_env[[1]]
    package_env <- unserialize(serialized_env)
    # Convert to an environment
    package_env <- list2env(package_env, parent = emptyenv())

    # Attach the environment to the search path
    attach(package_env, name = paste0("package:", package),
    warn.conflicts = FALSE
    )

    # # Attach the deserialized environment directly to `.GlobalEnv`
    # parent.env(package_env) <- .GlobalEnv
    # environment() <- package_env # cannot do this
    message(paste("Package", package, "successfully loaded from cache."))

  } else {
    stop(paste("No cached environment found for package:", package))
  }
}


.xx_load_environment_from_db <- function(db_path = "cache.sqlite") {
  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = db_path)

  # Query the most recent environment
  result <- DBI::dbGetQuery(con, "
    SELECT serialized_env FROM environment_cache
    ORDER BY timestamp DESC LIMIT 1
  ")

  DBI::dbDisconnect(con)

  if (nrow(result) > 0) {
    # Deserialize and load the environment
    loaded_env <- unserialize(result$serialized_env[[1]])
    # list2env(loaded_env, envir = .GlobalEnv)
    for (name in ls(loaded_env)) {
      assign(name, get(name, envir = loaded_env), envir = .GlobalEnv)
    }
    message("Environment successfully loaded from cache.")
  } else {
    stop("No cached environment found.")
  }
}


check_version_match_db <- function(db_path, packages_to_load, con = NULL, close = is.null(con)) {
  if (is.null(con)) {
    con <- DBI::dbConnect(RSQLite::SQLite(), dbname = db_path)
  }

  # Query all stored package versions
  stored_versions <- DBI::dbGetQuery(con, "
    SELECT package, version FROM package_cache
  ")
  if (close) {
    DBI::dbDisconnect(con)
    message("Connection closed.")
  }

  # Get current versions of requested packages
  installed_versions <- data.frame(
    package = packages_to_load,
    version = sapply(packages_to_load, function(pkg) as.character(packageVersion(pkg))),
    stringsAsFactors = FALSE
  )

  # Check if all requested packages are cached and versions match
  all_match <- all(installed_versions$package %in% stored_versions$package) &&
    all(installed_versions$version == stored_versions$version[match(installed_versions$package, stored_versions$package)])

  return(all_match)
}

# Check if versions in DB match planned packages
.xxcheck_version_match_with_packages <- function(db_path, packages_to_load) {
  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = db_path)
  stored_versions <- DBI::dbGetQuery(con, "SELECT * FROM package_cache")
  DBI::dbDisconnect(con)

  if (nrow(stored_versions) == 0) {
    return(FALSE)  # No versions in DB
  }

  # installed_versions <- sapply(packages_to_load, function(pkg) {
  #   as.character(packageVersion(pkg))
  # })

  # Get current versions of requested packages
  installed_versions <- data.frame(
    package = packages_to_load,
    version = sapply(packages_to_load, function(pkg) as.character(packageVersion(pkg))),
    stringsAsFactors = FALSE
  )
  # Check if all requested packages are cached and versions match
  all_match <- all(installed_versions$package %in% stored_versions$package) &&
    all(installed_versions$version == stored_versions$version[match(installed_versions$package, stored_versions$package)])

}

# Load packages and log versions to DB
._xxload_packages_and_log <- function(db_path, packages_to_load) {
  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = db_path)

  for (pkg in packages_to_load) {
    message(paste("Loading package:", pkg))
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }

  # Log package versions to DB
  installed_versions <- data.frame(
    package = packages_to_load,
    version = sapply(packages_to_load, function(pkg) {
      as.character(packageVersion(pkg))
    }),
    stringsAsFactors = FALSE
  )

  DBI::dbWriteTable(
    con, "package_versions", installed_versions,
    overwrite = TRUE, append = FALSE
  )

  # Serialize and save environment
  serialized_env <- serialize(environment(), NULL)
  serialized_env_blob <- I(list(serialized_env))

  version_hash <- digest::digest(installed_versions$version)
  DBI::dbExecute(
    con,
    "INSERT INTO environment_cache (serialized_env, version_hash) VALUES (?, ?)",
    params = list(serialized_env_blob, version_hash)
  )

  DBI::dbDisconnect(con)
}


handle_package_cache <- function(package, db_path, con = NULL, force_close = FALSE) {
  # Check if the package exists in the cache

  close <- is.null(con) || force_close
  if (is.null(con)) {
    con <- DBI::dbConnect(RSQLite::SQLite(), dbname = db_path)
  }

  result <- DBI::dbGetQuery(con, "
    SELECT version FROM package_cache
    WHERE package = ?
    LIMIT 1",
    params = list(package)
  )

  if (close) {
    DBI::dbDisconnect(con)
    message("Connection closed.")
    con <- NULL
  }


  # Get the current installed version
  current_version <- as.character(packageVersion(package))

  if (nrow(result) == 0) {
    # Package not in cache, load and store
    message(paste("Package", package, "not found in cache. Loading and storing."))
    save_package_to_db(package, db_path, con=con)
  } else if (result$version[[1]] != current_version) {
    # Version mismatch, reload and replace
    message(paste("Version mismatch for package", package, ". Reloading and replacing cache."))
    save_package_to_db(package, db_path, con=con)
  } else {
    # Version matches, load from cache
    message(paste("Version matches for package", package, ". Loading from cache."))
    load_package_from_db(package, db_path, con=con)
  }
}

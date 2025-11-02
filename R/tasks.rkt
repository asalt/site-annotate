;; experimental. not for use

#lang racket/base


(require racket/system          ; For `system*`
         racket/file            ; For `path`
         racket/path)           ; For `simplify-path`

;; Define a function to create a temporary file and execute it
(define (run-r-code-with-params params)
  (define r-folder (simplify-path (build-path (current-directory) ".." "R")))
  (define run-source (build-path r-folder "run.R"))

  ;; Check if the source file exists
  (unless (file-exists? run-source)
    (error 'run-r-code-with-params (format "File not found: ~a" run-source)))

  ;; Build the R code as a string
  (define r-code
    (string-join
     (append
      ;; Add parameters as R variables
      (for/list ([k (hash-keys params)]
                 [v (in-hash-values params)])
        (format "~a <- '~a'" k v))
      ;; Add the R logic
      (list
       (format "source('~a')" (path->string run-source))
       "# Interactive debugging support
if (exists('debug_mode') && debug_mode) {
    options(error = recover)
}

# Run the main function
print(paste0('Output dir is: ', output_dir))
run(
    data_dir = data_dir,
    output_dir = output_dir,
    config_file = config_file,
    gct_file = gct_file,
    save_env = save_env
)"))
     "\n"))

  ;; Create a temporary file for the R script
  (define temp-file (make-temporary-file "temp-r-script-" #:suffix ".R"))

  ;; Write the R code to the file
  (call-with-output-file temp-file
    (lambda (out)
      (fprintf out "~a" r-code))
    #:exists 'truncate)

  ;; Execute the R script
  (with-handlers ([exn:fail:system? (lambda (e)
                                      (displayln (exn-message e))
                                      (error 'run-r-code-with-params "Error running R script"))])
    (system* "Rscript" (path->string temp-file)))

  ;; Optionally delete the temporary file
  (delete-file temp-file))


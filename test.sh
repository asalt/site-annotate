test_py () {
    echo "running python tests"
    pytest ./tests/
}

test_R () {
    # Running R tests
    echo "Running R tests..."
    for f in ./tests/unit_tests/test_*R; do
        echo "Running $f..."
        (
            cd $(dirname $f)       # Switch to the directory where the test file is located
            Rscript $(basename $f) # Run the R script
        )
    done
}


test_R

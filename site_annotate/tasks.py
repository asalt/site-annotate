# tasks
import os
import subprocess
import tempfile
import pathlib
from rpy2 import robjects


def run_r_code_with_params(params_dict, interactive=False):
    # Generate R code dynamically
    r_code = "\n".join([f"{key} <- '{value}'" for key, value in params_dict.items()])
    run_source = os.path.join(pathlib.Path(__file__).parent, "..", "R", "run.R")

    run_source = pathlib.Path(__file__).parent / ".." / "R" / "run.R"
    r_folder = run_source.parent

    assert os.path.exists(run_source), f"File not found: {run_source}"

    # Escape curly braces in the R code
    r_code += f"""
    # Interactive debugging support
    options(error=recover)
    # if (exists("debug_mode") && debug_mode) {{
    #     options(error = recover)
    # }}
    setwd("{r_folder}")
    source("{run_source}")

    # Run the main function
    print(paste0('Output dir is: ', output_dir))
    run(
        data_dir = data_dir,
        output_dir = output_dir,
        config_file = config_file,
        gct_file = gct_file,
        save_env = save_env
    )
    """
    print(r_code)  # Optional: Print the R code for debugging

    # Create a temporary R script
    with tempfile.NamedTemporaryFile(
        suffix=".R", delete=False, mode="w"
    ) as temp_r_script:
        temp_r_script.write(r_code)
        temp_r_path = temp_r_script.name

    # Run the R script
    try:
        if not interactive:

            process = subprocess.Popen(
                ["R", "--no-save"],
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,  # Ensure text mode for stdin/stdout
            )
            # Send commands to the R session
            stdout, stderr = process.communicate(input=temp_r_path)
            subprocess.run(["R", "--no-save", "-f", temp_r_path], check=True)
            subprocess.run(["R", "--no-save", "-f", temp_r_path], check=True)
        else:
            robjects.r(r_code)

    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running R script: {e}")
    finally:
        # Optionally clean up the temporary file
        pathlib.Path(temp_r_path).unlink()


if __name__ == "__main__":
    # Example: Parameters to pass to the R script
    params_dict = {
        "data_dir": ".",
        "output_dir": "./output",
        "config_file": "config.yaml",
        "gct_file": "data.gct",
        "save_env": "FALSE",
        "debug_mode": "TRUE",  # Optional for debugging
    }

    run_r_code_with_params(params_dict)

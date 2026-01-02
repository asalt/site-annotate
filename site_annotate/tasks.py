# tasks
import os
import subprocess
import tempfile
import pathlib
import shutil
import textwrap
from typing import Optional, Dict, Any


def run_r_code_with_params(params_dict, interactive=False):
    def _r_literal(value: Any) -> str:
        if value is None:
            return "NULL"
        if isinstance(value, bool):
            return "TRUE" if value else "FALSE"
        if isinstance(value, (int, float)):
            return str(value)
        text = str(value)
        text = text.replace("\\", "\\\\").replace("'", "\\'")
        return f"'{text}'"

    run_source = pathlib.Path(__file__).parent / ".." / "R" / "run.R"
    run_source = run_source.resolve()
    r_folder = run_source.parent

    for key in ("data_dir", "config_file"):
        if key not in params_dict:
            raise ValueError(f"{key} not found")

    # Always define common parameters so the R call can safely reference them.
    params_dict = dict(params_dict)
    params_dict.setdefault("output_dir", None)
    params_dict.setdefault("gct_file", None)
    params_dict.setdefault("save_env", False)

    assert run_source.exists(), f"File not found: {run_source}"

    r_assignments = "\n".join([f"{key} <- {_r_literal(value)}" for key, value in params_dict.items()])
    r_code = f"""{r_assignments}

    # Interactive debugging support
    options(error=recover)
    setwd({_r_literal(str(r_folder))})
    source({_r_literal(str(run_source))})

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

    if interactive or os.environ.get("SITE_ANNOTATE_DEBUG_R"):
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
            subprocess.run(["Rscript", "--no-save", temp_r_path], check=True)
        else:
            try:
                from rpy2 import robjects
            except ImportError as exc:
                raise RuntimeError(
                    "Interactive mode requires rpy2; install it or disable --interactive"
                ) from exc

            robjects.r(r_code)

    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running R script: {e}")
    finally:
        # Optionally clean up the temporary file
        pathlib.Path(temp_r_path).unlink()


def maybe_compile_latex(summary_text: str, output_dir: str, basename: str = "summary", title: str = "Summary") -> str | None:
    """
    Optionally compile a simple LaTeX report using lualatex if available.
    Returns the path to the generated PDF or None if latex is unavailable or fails.
    """
    engine = shutil.which("lualatex") or shutil.which("pdflatex")
    if engine is None:
        return None

    outdir = pathlib.Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    tex_path = outdir / f"{basename}.tex"
    pdf_path = outdir / f"{basename}.pdf"

    # Try Jinja2 templating first, fallback to simple verbatim
    try:
        from jinja2 import Template  # type: ignore
        template_str = textwrap.dedent(r"""
            \documentclass[11pt]{article}
            \usepackage[margin=1in]{geometry}
            \usepackage{hyperref}
            \usepackage{lmodern}
            \title{ {{ title }} }
            \begin{document}
            \maketitle
            \begin{verbatim}
            {{ summary }}
            \end{verbatim}
            \end{document}
        """)
        tex = Template(template_str).render(title=title, summary=summary_text)
    except Exception:
        tex_tpl = textwrap.dedent(r"""
            \documentclass[11pt]{article}
            \usepackage[margin=1in]{geometry}
            \usepackage{hyperref}
            \usepackage{lmodern}
            \title{ {{TITLE}} }
            \begin{document}
            \maketitle
            \begin{verbatim}
            {{SUMMARY}}
            \end{verbatim}
            \end{document}
        """)
        tex = tex_tpl.replace("{{TITLE}}", title).replace("{{SUMMARY}}", summary_text)
    tex_path.write_text(tex)

    try:
        subprocess.run([engine, "-interaction=nonstopmode", "-halt-on-error", tex_path.name], cwd=str(outdir), check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        return None

    return str(pdf_path)


def compile_latex_from_template(context: Dict[str, Any], output_dir: str, basename: str, template: Optional[str] = None, title: str = "") -> Optional[str]:
    """Render a LaTeX template via Jinja2 if available and compile with lualatex/pdflatex."""
    engine = shutil.which("lualatex") or shutil.which("pdflatex")
    if engine is None:
        return None
    outdir = pathlib.Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    tex_path = outdir / f"{basename}.tex"
    pdf_path = outdir / f"{basename}.pdf"

    try:
        from jinja2 import Template  # type: ignore
        template_str = template or textwrap.dedent(r"""
            \documentclass[11pt]{article}
            \usepackage[margin=1in]{geometry}
            \usepackage{booktabs}
            \usepackage{lmodern}
            \title{ {{ title }} }
            \begin{document}
            \maketitle
            {% for section in sections %}
            \section*{ {{ section.title }} }
            \begin{tabular}{@{}{{ columns_spec }}@{}}
            \toprule
            {% for h in section.headers %} {{ h }} {% if not loop.last %}&{% endif %}{% endfor %} \\
            \midrule
            {% for row in section.rows %}
            {% for cell in row %} {{ cell }} {% if not loop.last %}&{% endif %}{% endfor %} \\
            {% endfor %}
            \bottomrule
            \end{tabular}
            {% endfor %}
            \end{document}
        """)
        tex = Template(template_str).render(**context)
    except Exception:
        return None

    tex_path.write_text(tex)
    try:
        subprocess.run([engine, "-interaction=nonstopmode", "-halt-on-error", tex_path.name], cwd=str(outdir), check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        return None
    return str(pdf_path)


def maybe_ollama_summarize(summary_text: str, output_dir: str, basename: str = "ollama_summary", model: str = "llama3.2:3b") -> Optional[str]:
    """If `ollama` is available, run a local model to summarize text and save output."""
    binary = shutil.which("ollama")
    if binary is None:
        return None
    outdir = pathlib.Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    prompt = (
        "Summarize the following results in 5-8 bullet points, "
        "highlighting significant hits and potential biological interpretations.\n\n" + summary_text
    )
    try:
        proc = subprocess.run([binary, "run", model], input=prompt.encode("utf-8"), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=True)
        out = proc.stdout.decode("utf-8", errors="replace")
    except subprocess.CalledProcessError:
        return None
    out_path = outdir / f"{basename}.md"
    out_path.write_text(out)
    return str(out_path)


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

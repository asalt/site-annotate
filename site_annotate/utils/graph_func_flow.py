import ast
import subprocess
from pathlib import Path
from collections import defaultdict
import glob

class CallGraphBuilder(ast.NodeVisitor):
    def __init__(self):
        self.graph = defaultdict(set)
        self.current_function = None
        self.constants = {}  # ✅ Fix: add this line


    def visit_FunctionDef(self, node):
        self.current_function = node.name
        self.generic_visit(node)
        self.current_function = None

    def visit_Call(self, node):
        if self.current_function is not None:
            if isinstance(node.func, ast.Name):
                callee = node.func.id
            elif isinstance(node.func, ast.Attribute):
                callee = node.func.attr
            else:
                callee = repr(node.func)

            self.graph[self.current_function].add(callee)
        self.generic_visit(node)


    def visit_Assign(self, node):
        if len(node.targets) != 1:
            return
        target = node.targets[0]
        if isinstance(target, ast.Name):
            var_name = target.id
            value = self._parse_value(node.value)
            self.constants[var_name] = value

    def _parse_value(self, node):
        if isinstance(node, ast.List):
            return [self._parse_value(elt) for elt in node.elts]
        elif isinstance(node, ast.Dict):
            return {
                self._parse_value(k): self._parse_value(v)
                for k, v in zip(node.keys, node.values)
            }
        elif isinstance(node, ast.Constant):
            return node.value
        elif isinstance(node, ast.Call):
            return f"CALL: {ast.unparse(node)}"
        else:
            return f"UNSUPPORTED: {ast.dump(node)}"


def build_call_graph(filename):
    with open(filename, "r") as f:
        tree = ast.parse(f.read())
    builder = CallGraphBuilder()
    builder.visit(tree)
    return builder.graph, builder.constants

def save_dot(graph, out_path: Path):
    all_funcs = set(graph.keys())
    all_callees = {callee for callees in graph.values() for callee in callees}
    all_nodes = all_funcs | all_callees

    def node_color(name):
        is_caller = name in graph
        is_callee = name in all_callees
        if is_caller and is_callee:
            return "gold"
        elif is_caller:
            return "lightblue"
        elif is_callee:
            return "lightgreen"
        else:
            return "gray"

    with open(out_path, "w") as f:
        f.write("digraph G {\n")
        f.write("    rankdir=LR;\n")
        f.write("    node [shape=box, style=filled];\n")

        # Write node styles
        for node in all_nodes:
            color = node_color(node)
            f.write(f'    "{node}" [fillcolor={color}];\n')

        # Write edges
        for caller, callees in graph.items():
            for callee in callees:
                f.write(f'    "{caller}" -> "{callee}" [arrowhead=vee, penwidth=1.5];\n')

        f.write("}\n")
    print(f"[✓] DOT saved to {out_path}")

def save_dot(graph, constants=None, out_path: Path = None):

    if constants is None:
        constants = dict()
    all_funcs = set(graph.keys())
    all_callees = {callee for callees in graph.values() for callee in callees}
    all_nodes = all_funcs | all_callees | set(constants.keys())

    def node_color(name):
        if name in constants:
            return "gray"
        is_caller = name in graph
        is_callee = name in all_callees
        if is_caller and is_callee:
            return "gold"
        elif is_caller:
            return "lightblue"
        elif is_callee:
            return "lightgreen"
        else:
            return "white"

    with open(out_path, "w") as f:
        f.write("digraph G {\n")
        f.write("    rankdir=LR;\n")
        f.write("    node [shape=box, style=filled];\n")

        for node in all_nodes:
            label = node
            if node in constants:
                label = f"{node}\\n[const]"
            f.write(f'    "{node}" [label="{label}", fillcolor={node_color(node)}];\n')

        for caller, callees in graph.items():
            for callee in callees:
                f.write(f'    "{caller}" -> "{callee}" [arrowhead=vee, penwidth=1.5];\n')

        f.write("}\n")
    print(f"[✓] DOT saved to {out_path}")

def render_dot_to_png_and_svg(dot_path: Path, out_png: Path):
    cmd = ["dot", "-Tpng", str(dot_path), "-o", str(out_png)]
    subprocess.run(cmd, check=True)
    print(f"[✓] PNG rendered to {out_png}")

    cmd = ["dot", "-Tsvg", str(dot_path), "-o", str(out_png).replace('png', 'svg')]
    subprocess.run(cmd, check=True)
    print(f"[✓] SVG rendered to {out_png.replace('png', 'svg')}")


def main():
    for source in glob.glob('../*py'):
        base = Path(source).stem
        graph, constants = build_call_graph(source)
        print(f"\n{'=' * 30} {base} {'=' * 30}")
        save_dot(graph, constants=None, out_path=f"{base}.dot")
        save_dot(graph, constants=constants, out_path=f"{base}_wconsts.dot")
        render_dot_to_png_and_svg(f"{base}_wconsts.dot", f"{base}_wconsts.png")

if __name__ == "__main__":
    main()


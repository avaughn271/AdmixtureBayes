"""running the premd script..."""

import sys
from call_AdmixtureBayes import main as run_main
from downstream_analysis_parser import run_posterior_main
from admixture_trees_to_consensus_tree import main as plot_main
from posterior_to_qpgraphs import main as qp_main
import argparse


def _collect_commands():
    commands = {}
    for name, value in list(globals().items()):
        if callable(value) and name.endswith("_command"):
            command_name = name.rsplit("_command", 1)[0]
            commands[command_name] = value
    return commands


# Commands
# For formatting main argument parseing message
class MixedFormatter(argparse.ArgumentDefaultsHelpFormatter,
                     argparse.RawDescriptionHelpFormatter):
    pass

def run_command(args):
    run_main(args)

def posterior_command(args):
    run_posterior_main(args)

def plot_command(args):
    plot_main(args)

def qp_command(args):
    qp_main(args)


# Main app
def main():
    "Main entry point for the script"
    commands = _collect_commands()
    commands_doc = "commands:\n{commands}".format(
        commands="\n".join(
            "  {name:10}:\t{doc}".format(name=name, doc=command.__doc__)
            for name, command in list(commands.items())
        )
    )

    parser = argparse.ArgumentParser(
        formatter_class=MixedFormatter,
        description="AdmixtureBayes program and auxillary tools",
        epilog=commands_doc
    )
    parser.add_argument(
        'command',
        help="Subcommand to run",
        metavar='command',
        choices=list(commands.keys())
    )
    parser.add_argument(
        'args', nargs="*",
        help="Arguments to subcommand"
    )

    args = parser.parse_args(sys.argv[1:2])  # only first two
    if args.command not in commands:
        parser.print_help()
        assert False,  "Unknown subcommand: {}".format(args.command)
        sys.exit(1)
    commands[args.command](sys.argv[2:])

# Run this in case the module is called as a program...
if __name__ == "__main__":
    main()

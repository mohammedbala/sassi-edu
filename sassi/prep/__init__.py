"""SASSI-EDU command interpreter ("UI" of ACS SASSI): ``.pre`` language, model editing, WRITE.

Quick start::

    from sassi.prep import Interpreter
    ui = Interpreter()
    ui.execute("N,1,0,0,0")
    ui.run_file("model.pre")          # INP
    ui.model                          # active sassi.model.SSIModel
    ui.sink.messages                  # messages (six classes, sassi.prep.messages.Kind)

See ``sassi/prep/README.md`` for adding commands.
"""
from .interpreter import Call, Interpreter, RunSummary, Variable, read_command_lines
from .lexer import LexError, parse_id_list, split_args, split_head
from .messages import Kind, Message, MessageSink, console_sink
from .registry import (CATALOGUE, CommandError, CommandReported, CommandSpec, all_commands, command, lookup,
                       register)

__all__ = ["Interpreter", "Call", "RunSummary", "Variable", "Kind", "Message", "MessageSink", "console_sink",
           "CommandError", "CommandReported", "CommandSpec", "command", "register", "lookup", "all_commands",
           "CATALOGUE", "LexError", "parse_id_list", "split_args", "split_head", "read_command_lines"]

##Wrap functions
def wrap_function(func, func_name):
    def class_initialization_wrap():
        def new_function(**kwargs):
            return_val= func(**kwargs)
            return {func_name:return_val},False
        return new_function
    return class_initialization_wrap


##wraps all functions ending in _custom_summary to return to downstream analysis parser for making custom functions.
def all_custom_summaries():
    commands = {}
    for name, value in list(globals().items()):
        if callable(value) and name.endswith("_custom_summary"):
            command_name = name.rsplit("_custom_summary", 1)[0]
            commands[command_name] = wrap_function(value, command_name)
    return commands



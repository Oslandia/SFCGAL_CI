import os

def define_env(env):
    """Define variables, macros, and filters.

    Args:
        env: The environment object used for defining macros and variables.

    Returns:
        None
    """
    @env.macro
    def include_file(filename: str) ->str:
        """Include the contents of a file.

        Args:
            filename (str): The name of the file to include.

        Returns:
            str: The content of the file with specific replacements made, 
                 or an error message if the file is not found.

        Raises:
            FileNotFoundError: If the specified file cannot be found.
        """
        try:
            with open(filename, "r", encoding="utf-8") as file:
                content = file.read().replace("website/docs/", "")
                content = content.replace("[AUTHORS](AUTHORS)", "[AUTHORS](./authors.md)")
                return content
        except FileNotFoundError:
            return f"File not found: {filename}"

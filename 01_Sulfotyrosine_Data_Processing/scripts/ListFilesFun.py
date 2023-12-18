import os

def getfiles(folder, extension, get_prefix=True, removal_strings=None):
    relative_paths = []
    prefixes = []
    # these are the filenames, so far i have come across these 3
    if removal_strings is None:
        removal_strings = ["_interact-ipro-ptm.pep.xml", "_interact-prob.pep.xml", "_interact-ipro.pep.xml"]

    for root, dirs, files in os.walk(folder):
        for file in files:
            if file.endswith(extension):
                relative_path = os.path.relpath(os.path.join(root, file), folder)
                relative_paths.append(relative_path)
                if get_prefix:
                    prefix = relative_path.replace("\\", "_")

                    # Remove specified strings if found in the prefix
                    for removal_string in removal_strings:
                        prefix = prefix.replace(removal_string, "")

                    # Strip the file extension if present
                    prefix, _ = os.path.splitext(prefix)

                    prefixes.append(prefix)

    return (relative_paths, prefixes) if get_prefix else relative_paths

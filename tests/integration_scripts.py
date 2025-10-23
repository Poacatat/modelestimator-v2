import re

def load_reference_output(output_string):

    output_string = output_string.splitlines(keepends=False)
    parsed_line =[]
    output_string = [line for line in output_string if line.strip() != '']
    for i in output_string:
        try:
            parts = [float(s) for s in re.split(r'\s+', i.strip())] 
        except ValueError:
            continue
        parsed_line.append(parts)

    return parsed_line

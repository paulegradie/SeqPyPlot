import os
from shutil import rmtree

def make_temp_upload_dir(directory):
    if os.path.exists(directory):
        rmtree(directory)
        os.mkdir(directory)
    else:
        os.mkdir(directory)


def retrieve_form_values(form, value_type):
    keys = [x for x in filter(lambda x: value_type in x, form)]
    values = [form.get(key) for key in keys]
    return keys, values


def safe_string_to_number(number):
    if '.' in number:
        try:
            number = float(number)
        except ValueError:
            pass
    else:
        try:
            number = int(number)
        except ValueError:
            pass
    return number

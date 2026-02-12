
# HTML Functions----------------------------------
def create_header(title, file_dict):
    html_string = f"""
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>{title}</title>
</head>
<body>
  <form>
    <select onchange="window.location.href = this.value;">
      <option value="">-- Select machine/config --</option>
"""
    for link, ff in file_dict.items():
        html_string += " "*6
        html_string += f'<option value="{ff}">{link}</option>\n'
    html_string += f"""
    </select>
  </form>
"""
    return html_string

def create_html_file(file_name, file_dict, title, content):
    html_string = create_header(title, file_dict)
    html_string += f"""
  <h1>{title}</h1>
"""
    for line in content:
        html_string += " "*4
        html_string += line+"\n"
    html_string += f"""
</body>
</html>
"""
    with open(file_name, 'w') as ff:
        ff.write(html_string)

def create_page_content(figs):
    content = []
    include_j = "cdn"
    for fig in figs:
        content.append(f"{fig.to_html(full_html=False, include_plotlyjs=include_j)}")
        include_j = False
    return content

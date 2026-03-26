import pydoc
import importlib.resources
import metapub.clinvarfetcher as module

css_data = importlib.resources.files('pydoc_data').joinpath('_pydoc.css').read_text()

pydoc.writedoc(module)  # this is the equivalent to pydoc -w
with open(module.__name__ + ".html") as inp:
    html = inp.read()
with open(module.__name__ + ".html", "w") as out:
    out.write(html.replace("</head>", "<style>%s</style></head>" % css_data))

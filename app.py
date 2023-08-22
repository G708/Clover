import dash
import pandas as pd
import pathlib
from dash import html
from dash import dcc

from dash.dependencies import Input, Output
from dash.exceptions import PreventUpdate

from apps.layout import html_layout


external_scripts = [
    {'src': 'https://ajax.googleapis.com/ajax/libs/jquery/1.11.2/jquery.min.js'},
    {'src': 'https://maxcdn.bootstrapcdn.com/bootstrap/3.2.0/js/bootstrap.min.js'}
]

app = dash.Dash(
    __name__,
    meta_tags=[{"name": "viewport", "content": "width=device-width, initial-scale=1"}],
)

app.title = 'Clover'

server = app.server
app.config.suppress_callback_exceptions = True
app.layout = html_layout()

# Path
BASE_PATH = pathlib.Path(__file__).parent.resolve()
DATA_PATH = BASE_PATH.joinpath("data").resolve()

# Read data
df = pd.read_csv(DATA_PATH.joinpath("DEPrior_gini_g2p.txt"), sep="\t", low_memory=False)

print(df.head())

# Run the server
if __name__ == "__main__":
    app.run_server(debug=True)
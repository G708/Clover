# Description: This file contains the components for the dashboard.

from dash import html
from dash import dcc
from datetime import datetime as dt

# Path
import pathlib
BASE_PATH = pathlib.Path(__file__).parent.resolve()
DATA_PATH = BASE_PATH.joinpath("data").resolve()


# control card have those options
# Tool description panel
# Controle panel
# 1. User file upload panel
# 2. Intaractive Ranking score setting panel (multiple scores)
# 3. submit / reset button
# 4. Intaractive 2D plot setting panel (X, Y, highlight)
# 5. Intaractive 3D plot setting panel (X, Y, Z, highlight)
# 6. Ranking fusion method setting panel
# 7. Ranking fusion rank setting panel (multiple ranks)
# 8. Word cloud setting panel
# 9. Tool version / link to github

# Tool description panel
def description_card():
    """
	Tool description panel

    :return: A Div containing dashboard title & descriptions.
    """
    return html.Div(
        id="description-card",
        children=[
            html.H5("Clover"),
            html.H3("Welcome to the Clover"),
            html.Div(
                id="intro",
                children="Explore DEG list by muptiple feature and sort by your own ranking score.",
            ),
        ],
    )


# Controle panel
def generate_control_card():
    """
	User Controle panel.
	There are 9 options in this panel.
	For each option, there is a title and a description.
	Discription and conponets are Collapse and Expand by clicking the title.

    :return: A Div containing controls for graphs.
    """
    return html.Div(
        id="control-card",
        children=[ 
			upload_panel()
        ],
    )

def upload_panel():
	"""
	User file upload panel
	
	:return: A Div containing controls for graphs.
	"""
	return html.Div(
		id="upload-panel",
		children=[
			html.H5("Upload your own data"),
			html.Div(
				id="upload",
				children=[
					dcc.Upload(
						id="upload-data",
						children=html.Div([
								'Drag and Drop or ',
								html.A('Select Files')
							]),
							style={
								'width': '100%',
								'height': '60px',
								'lineHeight': '60px',
								'borderWidth': '1px',
								'borderStyle': 'dashed',
								'borderRadius': '5px',
								'textAlign': 'center',
								'margin': '0px'
							},
						# Not allow multiple files to be uploaded
						multiple=False,
					),
				],
			),

		],
	)
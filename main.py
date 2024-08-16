from dash import Dash, html, dcc, Input, Output, State, callback_context
import dash_vtk
from dash_vtk.utils import to_mesh_state
import base64
import os
from web_stl_generator import generate_stl_from_params
from vtkmodules.vtkIOGeometry import vtkSTLReader
from vtkmodules.vtkFiltersSources import vtkPlaneSource

# Create a temporary directory to store uploaded files
UPLOAD_DIRECTORY = "uploads"
if not os.path.exists(UPLOAD_DIRECTORY):
    os.makedirs(UPLOAD_DIRECTORY)


# Function to clear the upload directory
def clear_upload_directory():
    for filename in os.listdir(UPLOAD_DIRECTORY):
        file_path = os.path.join(UPLOAD_DIRECTORY, filename)
        if os.path.isfile(file_path):
            os.remove(file_path)


# Dash setup
app = Dash(__name__, suppress_callback_exceptions=True)
server = app.server

app.layout = html.Div(
    style={
        "width": "100%",
        "height": "100%",
        "backgroundColor": "#f8f9fa",
        "fontFamily": "Arial, sans-serif",
        "padding": "20px"
    },
    children=[
        html.Div([
            html.H1("Crystal Lattice & Reciprocal Lattice Generator",
                    style={"textAlign": "center", "color": "#343a40"}),
            html.Div([
                html.Button("Reciprocal Lattice", id="reciprocal-lattice-btn", n_clicks=0,
                            style={"padding": "10px 20px", "margin": "5px", "backgroundColor": "#007bff",
                                   "color": "white", "border": "none", "borderRadius": "5px"}),
                html.Button("Crystal Lattice", id="crystal-lattice-btn", n_clicks=0,
                            style={"padding": "10px 20px", "margin": "5px", "backgroundColor": "#28a745",
                                   "color": "white", "border": "none", "borderRadius": "5px"}),
            ], style={"textAlign": "center", "margin": "20px 0"}),
        ], style={"boxShadow": "0 4px 8px 0 rgba(0,0,0,0.2)", "padding": "20px", "borderRadius": "10px",
                  "backgroundColor": "white", "margin": "20px auto", "width": "80%"}),
        html.Div(id="feature-content",
                 style={"margin": "10px auto", "padding": "20px", "backgroundColor": "white", "borderRadius": "10px",
                        "boxShadow": "0 4px 8px 0 rgba(0,0,0,0.2)", "width": "80%"})
    ],
)


@app.callback(
    Output("feature-content", "children"),
    [Input("reciprocal-lattice-btn", "n_clicks"),
     Input("crystal-lattice-btn", "n_clicks")]
)
def display_feature(rec_lattice_clicks, cry_lattice_clicks):
    if cry_lattice_clicks > rec_lattice_clicks:
        return [
            dcc.Upload(
                id="upload-cif",
                children=html.Div(["Drag and Drop or ", html.A("Select a CIF File")]),
                style={
                    "width": "100%",
                    "height": "60px",
                    "lineHeight": "60px",
                    "borderWidth": "1px",
                    "borderStyle": "dashed",
                    "borderRadius": "5px",
                    "textAlign": "center",
                    "margin": "10px 0",
                    "backgroundColor": "#e9ecef",
                },
                multiple=False,
            ),
            html.Div(id="output-upload", style={"margin": "10px 0", "color": "#495057"}),
            html.Div([
                html.Label(["Number of Unit Cells (x, y, z):", html.Span("?", title="Specify the number of unit cells in each direction (x, y, z).", style={"cursor": "pointer", "marginLeft": "5px", "color": "#007bff"})], style={"display": "block", "marginTop": "10px"}),
                dcc.Input(id="num-unit-cells-x", type="number", value=1,
                          style={"margin": "5px", "padding": "5px", "borderRadius": "5px",
                                 "border": "1px solid #ced4da", "width": "80px"}),
                dcc.Input(id="num-unit-cells-y", type="number", value=1,
                          style={"margin": "5px", "padding": "5px", "borderRadius": "5px",
                                 "border": "1px solid #ced4da", "width": "80px"}),
                dcc.Input(id="num-unit-cells-z", type="number", value=1,
                          style={"margin": "5px", "padding": "5px", "borderRadius": "5px",
                                 "border": "1px solid #ced4da", "width": "80px"}),
            ]),
            html.Div([
                html.Label(["Rotation Angles (x, y, z):", html.Span("?", title="Specify the rotation angles in degrees for the crystal in the x, y, and z directions.", style={"cursor": "pointer", "marginLeft": "5px", "color": "#007bff"})], style={"display": "block", "marginTop": "10px"}),
                dcc.Input(id="rotation-x", type="number", value=0,
                          style={"margin": "5px", "padding": "5px", "borderRadius": "5px",
                                 "border": "1px solid #ced4da", "width": "80px"}),
                dcc.Input(id="rotation-y", type="number", value=0,
                          style={"margin": "5px", "padding": "5px", "borderRadius": "5px",
                                 "border": "1px solid #ced4da", "width": "80px"}),
                dcc.Input(id="rotation-z", type="number", value=0,
                          style={"margin": "5px", "padding": "5px", "borderRadius": "5px",
                                 "border": "1px solid #ced4da", "width": "80px"}),
            ]),
            html.Div([
                html.Label(["Translation Vector (x, y, z):", html.Span("?", title="Specify the translation vector components for moving the crystal in the x, y, and z directions.", style={"cursor": "pointer", "marginLeft": "5px", "color": "#007bff"})], style={"display": "block", "marginTop": "10px"}),
                dcc.Input(id="translation-x", type="number", value=0,
                          style={"margin": "5px", "padding": "5px", "borderRadius": "5px",
                                 "border": "1px solid #ced4da", "width": "80px"}),
                dcc.Input(id="translation-y", type="number", value=0,
                          style={"margin": "5px", "padding": "5px", "borderRadius": "5px",
                                 "border": "1px solid #ced4da", "width": "80px"}),
                dcc.Input(id="translation-z", type="number", value=0,
                          style={"margin": "5px", "padding": "5px", "borderRadius": "5px",
                                 "border": "1px solid #ced4da", "width": "80px"}),
            ]),
            html.Div([
                html.Label(["Base Level:", html.Span("?", title="Specify the base level (height) for the crystal structure.", style={"cursor": "pointer", "marginLeft": "5px", "color": "#007bff"})], style={"display": "block", "marginTop": "10px"}),
                dcc.Input(id="base-level", type="number", value=0,
                          style={"margin": "5px", "padding": "5px", "borderRadius": "5px",
                                 "border": "1px solid #ced4da", "width": "80px"}),
            ]),
            html.Div([
                html.Label(["Is Primitive:", html.Span("?", title="Check this if you like to draw primitive lattice structure.", style={"cursor": "pointer", "marginLeft": "5px", "color": "#007bff"})], style={"display": "block", "marginTop": "10px"}),
                dcc.Checklist(id="is-primitive", options=[{'label': '', 'value': 'isPrimitive'}], value=[],
                              style={"margin": "5px", "padding": "5px"}),
            ]),
            html.Div([
                html.Label(["Target Atoms:", html.Span("?", title="Specify the target atoms, separated by commas.", style={"cursor": "pointer", "marginLeft": "5px", "color": "#007bff"})], style={"display": "block", "marginTop": "10px"}),
                dcc.Input(id="target-atoms", type="text", value=None,
                          style={"margin": "5px", "padding": "5px", "borderRadius": "5px",
                                 "border": "1px solid #ced4da", "width": "150px"}),
            ]),
            html.Div([
                html.Label(["Site Index Spin (format: index1:[x1,y1,z1],index2:[x2,y2,z2],...):", html.Span("?", title="Specify the spin configurations for site indices in the format provided.", style={"cursor": "pointer", "marginLeft": "5px", "color": "#007bff"})], style={"display": "block", "marginTop": "10px"}),
                dcc.Input(id="site-index-spin", type="text", value=None,
                          style={"margin": "5px", "padding": "5px", "borderRadius": "5px",
                                 "border": "1px solid #ced4da", "width": "300px"}),
            ]),
            html.Div([
                html.Label(["Tolerance:", html.Span("?", title="Specify the tolerance value for the STL generation.", style={"cursor": "pointer", "marginLeft": "5px", "color": "#007bff"})], style={"display": "block", "marginTop": "10px"}),
                dcc.Input(id="tolerance", type="number", value=0.1,
                          style={"margin": "5px", "padding": "5px", "borderRadius": "5px",
                                 "border": "1px solid #ced4da", "width": "80px"}),
            ]),
            html.Div([
                html.Label(["Add Supports:", html.Span("?", title="Check this if you want to add supports to the structure.", style={"cursor": "pointer", "marginLeft": "5px", "color": "#007bff"})], style={"display": "block", "marginTop": "10px"}),
                dcc.Checklist(id="add-supports-flag", options=[{'label': '', 'value': 'addSupports'}], value=[],
                              style={"margin": "5px", "padding": "5px"}),
            ]),
            html.Button("Generate STL", id="generate-stl", n_clicks=0,
                        style={"margin": "10px 0", "padding": "10px 20px", "backgroundColor": "#17a2b8",
                               "color": "white", "border": "none", "borderRadius": "5px"}),
            html.Button("Test Print", id="test-print-btn", n_clicks=0,
                        style={"margin": "10px 0", "padding": "10px 20px", "backgroundColor": "#ffc107",
                               "color": "black", "border": "none", "borderRadius": "5px"}),
            dcc.Download(id="download-stl"),
            html.Div(id="output-stl-path", style={"display": "none"}),
            html.Div(id="output-stl", style={"margin": "10px 0", "height": "400px"}),
            html.Button("Download STL", id="download-stl-btn", n_clicks=0,
                        style={"margin": "10px 0", "padding": "10px 20px", "backgroundColor": "#6c757d",
                               "color": "white", "border": "none", "borderRadius": "5px", "display": "none"}),
        ]
    elif rec_lattice_clicks > cry_lattice_clicks:
        return html.Div("Feature in progress", style={"color": "#6c757d", "textAlign": "center", "padding": "20px"})

    return "Select an option to proceed."

# The rest of your callbacks remain the same

if __name__ == "__main__":
    app.run_server(debug=True)

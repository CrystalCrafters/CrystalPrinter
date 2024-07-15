from dash import Dash, html, dcc, Input, Output, State
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
        "fontFamily": "Arial, sans-serif"
    },
    children=[
        html.Div([
            html.H1("Crystal Lattice & Reciprocal Lattice Generator", style={"textAlign": "center", "padding": "20px", "color": "#343a40"}),
            html.Div([
                html.Button("Reciprocal Lattice", id="reciprocal-lattice-btn", n_clicks=0, style={"padding": "10px 20px", "margin": "5px", "backgroundColor": "#007bff", "color": "white", "border": "none", "borderRadius": "5px"}),
                html.Button("Crystal Lattice", id="crystal-lattice-btn", n_clicks=0, style={"padding": "10px 20px", "margin": "5px", "backgroundColor": "#28a745", "color": "white", "border": "none", "borderRadius": "5px"}),
            ], style={"textAlign": "center", "margin": "20px 0"}),
        ], style={"boxShadow": "0 4px 8px 0 rgba(0,0,0,0.2)", "padding": "20px", "borderRadius": "10px", "backgroundColor": "white", "margin": "20px"}),
        html.Div(id="feature-content", style={"margin": "10px", "padding": "20px", "backgroundColor": "white", "borderRadius": "10px", "boxShadow": "0 4px 8px 0 rgba(0,0,0,0.2)"})
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
                    "margin": "10px",
                    "backgroundColor": "#e9ecef",
                },
                multiple=False,
            ),
            html.Div(id="output-upload", style={"margin": "10px", "color": "#495057"}),
            html.Div([
                html.Label("Number of Unit Cells (x, y, z):", style={"display": "block", "marginTop": "10px"}),
                dcc.Input(id="num-unit-cells-x", type="number", value=1, style={"margin": "5px", "padding": "5px", "borderRadius": "5px", "border": "1px solid #ced4da"}),
                dcc.Input(id="num-unit-cells-y", type="number", value=1, style={"margin": "5px", "padding": "5px", "borderRadius": "5px", "border": "1px solid #ced4da"}),
                dcc.Input(id="num-unit-cells-z", type="number", value=1, style={"margin": "5px", "padding": "5px", "borderRadius": "5px", "border": "1px solid #ced4da"}),
            ]),
            html.Div([
                html.Label("Rotation Angles (x, y, z):", style={"display": "block", "marginTop": "10px"}),
                dcc.Input(id="rotation-x", type="number", value=0, style={"margin": "5px", "padding": "5px", "borderRadius": "5px", "border": "1px solid #ced4da"}),
                dcc.Input(id="rotation-y", type="number", value=0, style={"margin": "5px", "padding": "5px", "borderRadius": "5px", "border": "1px solid #ced4da"}),
                dcc.Input(id="rotation-z", type="number", value=0, style={"margin": "5px", "padding": "5px", "borderRadius": "5px", "border": "1px solid #ced4da"}),
            ]),
            html.Div([
                html.Label("Translation Vector (x, y, z):", style={"display": "block", "marginTop": "10px"}),
                dcc.Input(id="translation-x", type="number", value=0, style={"margin": "5px", "padding": "5px", "borderRadius": "5px", "border": "1px solid #ced4da"}),
                dcc.Input(id="translation-y", type="number", value=0, style={"margin": "5px", "padding": "5px", "borderRadius": "5px", "border": "1px solid #ced4da"}),
                dcc.Input(id="translation-z", type="number", value=0, style={"margin": "5px", "padding": "5px", "borderRadius": "5px", "border": "1px solid #ced4da"}),
            ]),
            html.Div([
                html.Label("Base Level:", style={"display": "block", "marginTop": "10px"}),
                dcc.Input(id="base-level", type="number", value=0, style={"margin": "5px", "padding": "5px", "borderRadius": "5px", "border": "1px solid #ced4da"}),
            ]),
            html.Div([
                html.Label("Is Primitive:", style={"display": "block", "marginTop": "10px"}),
                dcc.Checklist(id="is-primitive", options=[{'label': '', 'value': 'isPrimitive'}], value=[], style={"margin": "5px", "padding": "5px"}),
            ]),
            html.Div([
                html.Label("Target Atoms:", style={"display": "block", "marginTop": "10px"}),
                dcc.Input(id="target-atoms", type="text", value=None, style={"margin": "5px", "padding": "5px", "borderRadius": "5px", "border": "1px solid #ced4da"}),
            ]),
            html.Div([
                html.Label("Site Index Spin:", style={"display": "block", "marginTop": "10px"}),
                dcc.Input(id="site-index-spin", type="text", value=None, style={"margin": "5px", "padding": "5px", "borderRadius": "5px", "border": "1px solid #ced4da"}),
            ]),
            html.Div([
                html.Label("Tolerance:", style={"display": "block", "marginTop": "10px"}),
                dcc.Input(id="tolerance", type="number", value=0.1, style={"margin": "5px", "padding": "5px", "borderRadius": "5px", "border": "1px solid #ced4da"}),
            ]),
            html.Div([
                html.Label("Add Supports:", style={"display": "block", "marginTop": "10px"}),
                dcc.Checklist(id="add-supports-flag", options=[{'label': '', 'value': 'addSupports'}], value=[], style={"margin": "5px", "padding": "5px"}),
            ]),
            html.Button("Generate STL", id="generate-stl", n_clicks=0, style={"margin": "10px", "padding": "10px 20px", "backgroundColor": "#17a2b8", "color": "white", "border": "none", "borderRadius": "5px"}),
            dcc.Download(id="download-stl"),
            html.Div(id="output-stl-path", style={"display": "none"}),
            html.Div(id="output-stl", style={"margin": "10px", "height": "400px"}),
            html.Button("Download STL", id="download-stl-btn", n_clicks=0, style={"margin": "10px", "padding": "10px 20px", "backgroundColor": "#6c757d", "color": "white", "border": "none", "borderRadius": "5px", "display": "none"}),
        ]
    elif rec_lattice_clicks > cry_lattice_clicks:
        return html.Div("Feature in progress", style={"color": "#6c757d", "textAlign": "center", "padding": "20px"})

    return "Select an option to proceed."

@app.callback(
    Output("output-upload", "children"),
    Input("upload-cif", "contents"),
    State("upload-cif", "filename"),
)
def save_upload(contents, filename):
    if contents is not None:
        clear_upload_directory()
        data = contents.encode("utf8").split(b";base64,")[1]
        cif_path = os.path.join(UPLOAD_DIRECTORY, filename)
        with open(cif_path, "wb") as fp:
            fp.write(base64.decodebytes(data))
        return f"Uploaded file: {filename}"
    return "No file uploaded yet."

@app.callback(
    [Output("output-stl-path", "children"),
     Output("download-stl-btn", "style")],
    Input("generate-stl", "n_clicks"),
    State("upload-cif", "filename"),
    State("num-unit-cells-x", "value"),
    State("num-unit-cells-y", "value"),
    State("num-unit-cells-z", "value"),
    State("rotation-x", "value"),
    State("rotation-y", "value"),
    State("rotation-z", "value"),
    State("translation-x", "value"),
    State("translation-y", "value"),
    State("translation-z", "value"),
    State("base-level", "value"),
    State("is-primitive", "value"),
    State("target-atoms", "value"),
    State("site-index-spin", "value"),
    State("tolerance", "value"),
    State("add-supports-flag", "value"))
def generate_stl(n_clicks, filename, num_x, num_y, num_z, rot_x, rot_y, rot_z, trans_x, trans_y, trans_z, base_level, is_primitive, target_atoms, site_index_spin, tolerance, add_supports_flag):
    if n_clicks > 0 and filename:
        cif_path = os.path.join(UPLOAD_DIRECTORY, filename)
        num_unit_cells = [num_x, num_y, num_z]
        rotation_angles = [rot_x, rot_y, rot_z]
        translation_vector = [trans_x, trans_y, trans_z]
        is_primitive = True if 'isPrimitive' in is_primitive else False
        target_atoms = target_atoms.split(',') if target_atoms else None
        site_index_spin = site_index_spin.split(',') if site_index_spin else None
        add_supports_flag = True if 'addSupports' in add_supports_flag else False

        stl_file_path = generate_stl_from_params(
            cif_path,
            num_unit_cells,
            rotation_angles,
            translation_vector,
            base_level,
            is_primitive,
            target_atoms,
            site_index_spin,
            tolerance,
            add_supports_flag
        )

        if os.path.exists(stl_file_path):
            return stl_file_path, {"display": "block"}
        else:
            return "Failed to generate STL file.", {"display": "none"}

    return "Click 'Generate STL' to create the STL file.", {"display": "none"}

@app.callback(
    Output("download-stl", "data"),
    Input("download-stl-btn", "n_clicks"),
    State("output-stl-path", "children"),
    prevent_initial_call=True
)
def download_stl(n_clicks, stl_file_path):
    if n_clicks > 0 and stl_file_path and os.path.exists(stl_file_path):
        return dcc.send_file(stl_file_path)
    return None

@app.callback(
    Output("output-stl", "children"),
    Input("output-stl-path", "children"),
    State("base-level", "value"),
)
def render_stl(stl_file_path, base_level):
    if stl_file_path and os.path.exists(stl_file_path):
        stl_reader = vtkSTLReader()
        stl_reader.SetFileName(stl_file_path)
        stl_reader.Update()

        dataset = stl_reader.GetOutput()

        if dataset is None:
            return "Failed to read STL file."

        mesh_state = to_mesh_state(dataset)

        # Compute the bounds of the STL model
        bounds = dataset.GetBounds()
        x_min, x_max, y_min, y_max, z_min, z_max = bounds
        x_size = x_max - x_min
        y_size = y_max - y_min

        # Create the plane with the computed size
        plane_source = vtkPlaneSource()
        plane_source.SetOrigin(x_min, y_min, base_level)
        plane_source.SetPoint1(x_max, y_min, base_level)
        plane_source.SetPoint2(x_min, y_max, base_level)
        plane_source.Update()

        plane_dataset = plane_source.GetOutput()
        plane_mesh_state = to_mesh_state(plane_dataset)

        content = dash_vtk.View([
            dash_vtk.GeometryRepresentation([
                dash_vtk.Mesh(state=plane_mesh_state),
            ], property={"color": [0.8, 0.8, 0.8]}),
            dash_vtk.GeometryRepresentation([
                dash_vtk.Mesh(state=mesh_state),
            ])
        ])

        return content

    return "No STL file to display."

if __name__ == "__main__":
    app.run_server(debug=True)

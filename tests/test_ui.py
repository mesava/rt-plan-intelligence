"""Smoke test of the real Streamlit app using a synthetic-only example."""
from pathlib import Path
from streamlit.testing.v1 import AppTest

def test_streamlit_opens_and_demo_renders():
    root=Path(__file__).resolve().parents[1]
    app=AppTest.from_file(str(root/"app.py"),default_timeout=30).run()
    assert len(app.exception)==0, [x.message for x in app.exception]
    assert len(app.checkbox)>0
    app.checkbox[0].set_value(True).run()
    assert len(app.exception)==0, [x.message for x in app.exception]
    assert len(app.dataframe)>0

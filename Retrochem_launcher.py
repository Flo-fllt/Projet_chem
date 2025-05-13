import sys
import os
import streamlit.web.cli as stcli

def main():
    script_path = os.path.join(os.path.dirname(__file__), "RetroChem", "Retrochem.py")
    sys.argv = ["streamlit", "run", script_path]
    sys.exit(stcli.main())

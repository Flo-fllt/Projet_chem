# retrochem_launcher.py

def main():
    import sys
    import streamlit.web.cli as stcli
    sys.argv = ["streamlit", "run", "Scripts/Retrochem.py"]
    sys.exit(stcli.main())

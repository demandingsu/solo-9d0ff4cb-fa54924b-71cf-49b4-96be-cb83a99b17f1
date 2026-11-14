import sys
from pathlib import Path

# Allow `python -m pytest` from anywhere in the repository.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

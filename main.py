"""
Created: 15 Oct 2022
Project: Cell-type marker heritability
"""
from CellMarkerHandler import CellTypeMarker


def main():
    cell_type_marker = CellTypeMarker(path_to_tsv='./Cell-type marker/PanglaoDB_markers_27_Mar_2020.tsv')


if __name__ == "__main__":
    main()

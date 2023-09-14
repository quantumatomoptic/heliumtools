import importlib, sys


def main():
    module_name = "figure0"  # Nom du module que vous souhaitez importer
    sys.path.append("/home/victor/heliumtools/heliumtools/GUIs/base")
    try:
        mon_module = importlib.import_module(module_name)
        mon_module.ma_fonction()  # Appeler une fonction spécifique du module
    except ImportError:
        print(f"Impossible d'importer le module {module_name}")


if __name__ == "__main__":
    main()

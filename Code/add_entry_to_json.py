import json
import sys

def add_entry_to_json(file_path, key, value):
    try:
        # Read the existing data from the JSON file
        with open(file_path, 'r') as file:
            data = json.load(file)
    except FileNotFoundError:
        # If the file does not exist, start with an empty dictionary
        data = {}
    except json.JSONDecodeError:
        # If the file is not a valid JSON, start with an empty dictionary
        data = {}

    # Add the new entry
    data[key] = value

    # Write the updated data back to the JSON file
    with open(file_path, 'w') as file:
        json.dump(data, file, indent=4)

    print(f"Added entry: {key} -> {value} to {file_path}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python add_entry_to_json.py <file_path> <key> <value>")
        sys.exit(1)

    file_path = sys.argv[1]
    key = sys.argv[2]
    value = sys.argv[3]

    add_entry_to_json(file_path, key, value)

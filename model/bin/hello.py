import sys

def main():
    if len(sys.argv) < 2:
        print("No argument provided!")
        return

    arg = sys.argv[1]
    print(f"Hello from Python! You passed the argument: {arg}")

if __name__ == "__main__":
    main()

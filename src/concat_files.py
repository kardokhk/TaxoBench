import sys

def main():
    for file in sys.argv[1:]:
        with open(file, "r") as f:
            for line in f:
                line = line.rstrip('\r\n')
                sys.stdout.write(line + '\n')

if __name__ == "__main__":
    main()

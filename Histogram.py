import matplotlib.pyplot as plt

def get_numbers(file):
     fin = open(file, "r")
     for line in fin:
          numbers.append(line.strip())
     print numbers

numbers = get_numbers(sys.argv[1])

plt.hist(numbers, bins = 40)
plt.xlabel("Value")
plt.ylabel("Frequency")
plt.show()


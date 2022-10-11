from decimal import Decimal
from sys import argv
from unicodedata import name

## The array stats class encapsulate the main parameters to make statics on the result
## and provide a method to print the results
class ArrayStats:
    # Constructor
    def __init__(self, number_wrong_values, total_values, percentage_wrong_values):
        self.number_wrong_values = number_wrong_values
        self.total_values = total_values
        self.percentage_wrong_values = percentage_wrong_values
    
    # Print the stats given in input the name of the array
    def print_stats(self, array_name):
        print("")
        print(f"Vector {array_name}")
        print(f"Number of values analyzed: {self.total_values}")
        print(f"Number of different values: {self.number_wrong_values}")
        print(f"Percentage of different values: {self.percentage_wrong_values}")

## The FileAnalyzer class allows reading data from file passed as argument
## and a method wich return a list of statics object that will be used
## to print out the result
class FileAnalyzer:
    def __read_arrays_from_file(self, file_name):
        # declare the output list
        output = []
        # open the file in read mode
        file = open(file_name, "r")
        # save the number of arrays in the file
        num_arrays = int(file.readline())
        # save each array in a temporary list and then
        # append to the final output list
        for i in range(0, num_arrays):
            n = int(file.readline())
            temp_list = []
            for i in range(0, n):
                value = Decimal(file.readline())
                temp_list.append(value)
            output.append(list(temp_list))
        # close the file
        file.close()
        return output

    def analize_output_files(self, file_name1, file_name2):
        results = []
        outputC = self.__read_arrays_from_file(file_name1)
        outputMATLAB = self.__read_arrays_from_file(file_name2)
        ## Checks
        # Check the number of arrays
        if (len(outputC) != len(outputMATLAB)):
            print("The array number of the two files is not the same")
            return
        # Check the number of element of corrisponding arrays
        for i in range(0, len(outputC)):
            if (len(outputC[i]) != len(outputMATLAB[i])):
                print("The array number of the %d-th array is not the same" % i)
                return
        # Built the result lsit
        for i in range(0, len(outputC)):
            # Build a list in which each value is a difference between the value from
            # the list outputC and the value from the list outputMATLAB
            differences = list(map(lambda x, y: y - x, outputC[i], outputMATLAB[i]))
            # Compute the number of wrong values, which means the values in difference
            # that are not equal to zero
            num_wrong_values = len(differences) - differences.count(0)
            # Calculate the percentage of the wrong values
            percentage_wrong_values = (num_wrong_values / len(differences)) * 100
            # Encapsulate the calculated data in an object and save it in a list
            local_stats = ArrayStats(num_wrong_values, len(differences), percentage_wrong_values)
            results.append(local_stats)
        return results

if __name__ == "__main__":
    # Create a new FileAnalyzer object
    file_analyzer = FileAnalyzer()
    # Call method to analyze each file
    total_stats = file_analyzer.analize_output_files(argv[1], argv[2])
    # Print results
    total_stats[0].print_stats("y")
    total_stats[1].print_stats("yT")
    total_stats[2].print_stats("t")
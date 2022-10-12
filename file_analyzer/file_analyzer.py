from collections import Counter
from decimal import Decimal
from sys import argv

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
            differences = list(map(lambda x, y: abs(y - x), outputC[i], outputMATLAB[i]))
            # Build a dictonary that for each value in the list associates the number
            # of occourrences of that value
            local_stats = dict(Counter(differences))
            # Append to the result a list in which the first element is the number of elements
            # and the second is the dictonary of occourrences
            results.append([len(differences), local_stats])
        return results

def print_local_stats(array_name, local_list):
    print("")
    print(f"Vector {array_name}")
    print(f"Number of values analyzed: {local_list[0]}")
    print("difference -> #values")
    for key in local_list[1]:
        print(f"{key} -> {local_list[1][key]} about {(local_list[1][key] / local_list[0] * 100):.2f}%")

if __name__ == "__main__":
    # Create a new FileAnalyzer object
    file_analyzer = FileAnalyzer()
    # Call method to analyze each file
    total_stats = file_analyzer.analize_output_files(argv[1], argv[2])
    # Print results
    print_local_stats("y", total_stats[0])
    print_local_stats("yT", total_stats[1])
    print_local_stats("t", total_stats[2])
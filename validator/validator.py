from decimal import Decimal
from sys import argv
from unicodedata import name

class Validator:
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
        outputC = self.__read_arrays_from_file(file_name1)
        outputMATLAB = self.__read_arrays_from_file(file_name2)
        ## Checks
        # Check the number of arrays
        if (len(outputC) != len(outputMATLAB)):
            print("The array number of the two files is not the same")
            return
        ## Check the number of element of corrisponding arrays
        for i in range(0, len(outputC)):
            if (len(outputC[i]) != len(outputMATLAB[i])):
                print("The array number of the %d-th array is not the same" % i)
                return
        

if __name__ == "__main__":
    validator = Validator()
    validator.analize_output_files(argv[1], argv[2])
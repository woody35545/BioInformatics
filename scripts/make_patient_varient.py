import random as r
import copy

SIZE_OF_PATIENT_DATA = 122 
SIZE_OF_VARIANT = 56 
SEQUENCE = ['A','C','G','T','D']
index = [0] * SIZE_OF_VARIANT


def init_index():
  global index
  for i in range(SIZE_OF_VARIANT):
    index[i] = r.randrange(0,5600)

def make_noise(_original_data, num_of_varient): 
  original_list = []
  for i in range(5600):
    original_list.append(_original_data[i])
  
  for i in range(num_of_varient):
      init_index()
      res = copy.copy(original_list)
      for i in range(SIZE_OF_VARIANT):
        res[int(index[i])] = r.choice(SEQUENCE)
      f = open("../data/patient_noise.txt", 'a')
      f.write(''.join(res))
      f.write("\n")
      f.close()
  return res


for j in range(SIZE_OF_PATIENT_DATA):
    patient_file_path = "../exon/exon_SRR7785" + str(j+619) + ".txt"
    patient_file = open(patient_file_path, 'r')
    original_data = patient_file.readline()
    patient_file.close()
    make_noise(original_data,10)


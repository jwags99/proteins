def two_character_count_map(sequences):
  count_map = {}
  for sublist in sequences: # One list per file
    for protein in sublist:
      for i in range(len(str(protein)) - 1):
        two_char_str = protein[i] + protein[i+1]
        if two_char_str in count_map:
          count_map[two_char_str] += 1
        else:
          count_map[two_char_str] = 1
  return count_map

def three_character_count_map(sequences):
    count_map = {}
    for sublist in sequences: # One list per file
        for protein in sublist:
            for i in range(len(str(protein)) - 2):
                three_char_str = protein[i] + protein[i+1] + protein[i+2]
                if three_char_str in count_map:
                    count_map[three_char_str] += 1
                else:
                    count_map[three_char_str] = 1
    return count_map

def four_character_count_map(sequences):
  count_map = {}
  for sublist in sequences: # One list per file
    for protein in sublist:
      for i in range(len(str(protein)) - 3):
        four_char_str = protein[i] + protein[i+1] + protein[i+2] + protein[i+3]
        if four_char_str in count_map:
          count_map[four_char_str] += 1
        else:
          count_map[four_char_str] = 1
  return count_map

def two_char(protein_sequence): #Can make an array of arrays w/ bi_char and get rid of twoletters
  twoletters = []
  for i, file in enumerate(protein_sequence):
    # print(i, c)
    if (i != len(protein_sequence) - 1):
      for protein in file:
        for index in range(len(protein)):
          if (index != len(protein)-1):  
            bi_char = protein[index] + protein[index + 1]
            twoletters.append(bi_char)
          else:
            continue
  return twoletters

def three_char(protein):  #Can make an array of arrays w/ tri_char and get rid of twoletters
  threeletters = []
  for i, c in enumerate(protein):
    if (i != len(protein) - 1):
      for x in c:
        for index in range(len(x)):
          if (index != len(x)-2):  
            tri_char = x[index] + x[index + 1] + x[index + 2]
            threeletters.append(tri_char)
          else:
            break
  return threeletters

def four_char(protein):  #Can make an array of arrays w/ bi_char and get rid of twoletters
  fourletters = []
  for i, c in enumerate(protein):
    if (i != len(protein) - 1):
      for x in c:
        for index in range(len(x)):
          if (index != len(x)-3):  
            charmander_char = x[index] + x[index + 1] + x[index + 2] + x[index + 3]
            fourletters.append(charmander_char)
          else:
            break
  return fourletters

def character_count_map(sequences):
  count_map = {}
  for sublist in sequences: # One list per file
    for protein in sublist:
      for c in protein:
        if c in count_map:
          count_map[c] += 1
        else:
          count_map[c] = 1

  return sort_dictionary(count_map)

def sort_dictionary(orig_dictionary):
  sorted_list = sorted(orig_dictionary, key=orig_dictionary.get, reverse=True)
  sorted_dictionary = {}
  
  for i in sorted_list:
    sorted_dictionary[i] = orig_dictionary[i]

  return sorted_dictionary
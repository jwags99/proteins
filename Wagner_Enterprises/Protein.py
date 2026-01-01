class Protein:
    def __init__(self, sequence, name):
        self.seq = sequence
        self.name = name
        self.single_character_dict = {}
        self.two_character_dict = {}
        self.three_character_dict = {}
        self.four_character_dict = {}
        self.single_character_count_map()
        self.two_character_count_map()
        self.three_character_count_map()
        self.four_character_count_map()

    def single_character_count_map(self):
        for char in self.seq:
            if char in self.single_character_dict:
                self.single_character_dict[char] += 1
            else:
                self.single_character_dict[char] = 1
    
    def two_character_count_map(self):
        sequence_str = ''.join(self.seq)  # Join all lines/strings in the list into one continuous string
        for i in range(len(sequence_str) - 1):
            two_char_str = sequence_str[i] + sequence_str[i+1]
            if two_char_str in self.two_character_dict:
                self.two_character_dict[two_char_str] += 1
            else:
                self.two_character_dict[two_char_str] = 1
        # for i in range(len(str(self.seq)) - 1):
        #     two_char_str = self.seq[i] + self.seq[i+1]
        #     if two_char_str in self.two_character_dict:
        #         self.two_character_dict[two_char_str] += 1
        #     else:
        #         self.two_character_dict[two_char_str] = 1

    def three_character_count_map(self):
        sequence_str = ''.join(self.seq)  # Join all lines/strings in the list into one continuous string
        for i in range(len(sequence_str) - 2):
            three_char_str = sequence_str[i] + sequence_str[i+1] + sequence_str[i+2]
            if three_char_str in self.three_character_dict:
                self.three_character_dict[three_char_str] += 1
            else:
                self.three_character_dict[three_char_str] = 1
        # for i in range(len(str(self.seq)) - 2):
        #     three_char_str = self.seq[i] + self.seq[i+1] + self.seq[i+2]
        #     if three_char_str in self.three_character_dict:
        #         self.three_character_dict[three_char_str] += 1
        #     else:
        #         self.three_character_dict[three_char_str] = 1

    def four_character_count_map(self):
        sequence_str = ''.join(self.seq)  # Join all lines/strings in the list into one continuous string
        for i in range(len(sequence_str) - 3):
            four_char_str = sequence_str[i] + sequence_str[i+1] + sequence_str[i+2] + sequence_str[i+3]
            if four_char_str in self.four_character_dict:
                self.four_character_dict[four_char_str] += 1
            else:
                self.four_character_dict[four_char_str] = 1
        # for i in range(len(str(self.seq)) - 3):
        #     four_char_str = self.seq[i] + self.seq[i+1] + self.seq[i+2] + self.seq[i+3]
        #     if four_char_str in self.four_character_dict:
        #         self.four_character_dict[four_char_str] += 1
        #     else:
        #         self.four_character_dict[four_char_str] = 1
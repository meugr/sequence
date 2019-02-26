# ver 1.1
# Получение последовательности нуклеотидов между двумя заданными
# последовательностями. Обрабатывает заданную, комплиментарную и
# обратные им последовательности.

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data import IUPACData
from os import listdir, path


def get_sequence(path_file):
    """
    Получаем последовательность нуклеотидов из файла.
    Возвращаем строку содержащую последовательность, например 'AGAGAAAGATAG'
    """
    with open(path_file, "rb") as input_file:
        for i in SeqIO.parse(input_file, "abi"):
            sequence = Seq("".join(i[:]),)
        return sequence


def make_complimentary(sequence):
    """
    Создает комплиментарную последовательность принятой последовательности.
    Возвращает строку содержащую последовательность, например 'AGAGAAAGATAG'
    Комплиментарные нуклеотиды: a-t c-g r-y k-m n-n x-x b-v d-h
    """
    return sequence.complement()


def search(seq, subseq):
    """Находит в последовательности с неопределенностями нужную
    фланкирующую последовательность, возвращает ее и порядковый номер ее
    начала, начиная с нуля.
    """
    dataset = []
    subseq = list(subseq)
    for nt in seq:
        value = IUPACData.ambiguous_dna_values[nt]
        dataset.append(set(value))

    for start in range(len(dataset)):
        for offset in range(len(subseq)):
            try:
                if subseq[offset] not in dataset[start + offset]:
                    break
            except IndexError:
                break
        else:
            return start


def find_sequence(sequence, sequence_reverse,
                  comp_sequence, comp_sequence_reverse,
                  start, end):
    sequences = (sequence, sequence_reverse, comp_sequence,
                 comp_sequence_reverse)
    seq_type = -1
    for i in sequences:
        seq_type += 1

        seq_start = search(str(i), start)
        seq_end = search(str(i), end)
        if seq_start and seq_end:
            total_cou = 0
            for c in ('A', 'C', 'T', 'G'):
                cou = i[seq_start: seq_end + len(end)].count(c)
                total_cou += cou
            final_seq = i[seq_start: seq_end + len(end)]
            return final_seq, seq_start, seq_end,\
                seq_type, len(final_seq) - total_cou


def write_in_output_file(output_sequence_list):
    """
    Запись списка полученых в find_sequence() последовательностей в файл
    output.fasta
    """
    with open("output.fasta", "w") as file:
        SeqIO.write(output_sequence_list, file, "fasta")


def write_in_error_file(output_error_list):
    """
    Запись имени файла в котором произошла ошибка в файл errors.txt
    Возможные ошибки:
    1) отсутствие, искажение или дублирование последовательностей start и end
    2) Файл поврежден и не открывается
    """
    with open("error.txt", "w") as file:
        for i in output_error_list:
            file.write(i + '\n')


# +++++++++++++++++++++++++++++++MAIN+++++++++++++++++++++++++++++++

print("""
Получение последовательности нуклеотидов между двумя заданными
последовательностями. Обрабатывает заданную, комплиментарную и
обратные им последовательности. Учитывает неоднозначные нуклеотиды.

КАК ПОЛЬЗОВАТЬСЯ
Создайте в директории с программой директорию 'input', куда положите
Ваши файлы с расширением *.ab1
Введите начальную и конечную последовательность
Результат будет в файле output.fasta, список необработанных файлов
в error.txt

id последовательности - имя файла.
description - длина вставки, положение в плазмиде,
прямая, реверснутая, комплиментарная или комплиментарнвя реверснутая,
количество неопределенностей
""")

subdir = 'input'
start = input('Введите начало последовательности: ').upper()
end = input('Введите конец последовательности: ').upper()
output_sequence_list = []
output_error_list = []
filename_list = listdir(subdir)
seq_types = ('sequence', 'sequence_reverse', 'comp_sequence',
             'comp_sequence_reverse')

for filename in filename_list:
    path_file = path.join(subdir, filename)
    # Получаем 4 вида одной последовательности
    sequence = get_sequence(path_file)
    sequence_reverse = sequence[::-1]
    comp_sequence = make_complimentary(sequence)
    comp_sequence_reverse = comp_sequence[::-1]

    sequence_find_info = find_sequence(sequence, sequence_reverse,
                                       comp_sequence, comp_sequence_reverse,
                                       start, end)

    if sequence_find_info:
        # сбор данных для description
        length = len(sequence_find_info[0])
        ambiguity = sequence_find_info[4]
        description = f'length: {length}, \
position in the plasmid: {sequence_find_info[1:3]}, seq_type: \
{seq_types[sequence_find_info[3]]}, ambiguity: {ambiguity}'
        # Записываем последовательность в список.
        output_sequence_list.append(SeqRecord(sequence_find_info[0],
                                              id=filename,
                                              description=description))
        print(sequence_find_info[0])
    else:
        # Записываем название файла, в котором возникла ошибка.
        output_error_list.append(filename)

write_in_output_file(output_sequence_list)
write_in_error_file(output_error_list)
print(f'''
=========================================================
Готово!
Обработано файлов: {len(filename_list)}
Получено последовательностей: {len(output_sequence_list)}
''')

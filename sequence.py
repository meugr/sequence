# Получение последовательности нуклеотидов между двумя заданными
# последовательностями. Обрабатывает заданную, комплиментарную и
# обратные им последовательности.

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from os import listdir, path


def get_sequence(path_file):
    """
    Получаем последовательность нуклеотидов из файла.
    Возвращаем строку содержащую последовательность, например 'AGAGAAAGATAG'
    """
    with open(path_file, "rb") as input_file:
        for i in SeqIO.parse(input_file, "abi"):
            sequence = Seq("".join(i[:]))
        return sequence


def make_complimentary(sequence):
    """
    Создает комплиментарную последовательность принятой последовательности.
    Возвращает строку содержащую последовательность, например 'AGAGAAAGATAG'
    Комплиментарные нуклеотиды: a-t c-g r-y k-m n-n x-x b-v d-h
    """
    return sequence.complement()


def find_sequence(sequence, sequence_reverse,
                  comp_sequence, comp_sequence_reverse,
                  start, end):
    """
    Находит последовательность, заключенную между двумя данными - start и end.
    Возвращает строку-последовательность нуклеотидов
    """
    sequences = [sequence, sequence_reverse, comp_sequence,
                 comp_sequence_reverse]
    for i in sequences:
        if start in i and end in i:
            if i.count(start) == 1 and i.count(end) == 1:
                return (i.split(start))[1].split(end)[0]


def write_in_output_file(output_sequence_list):
    """
    Запись списка полученых в find_sequence() последовательностей в файл
    output.fasta
    """
    with open("output.fasta", "w") as file:
        SeqIO.write(output_sequence_list, file, "fasta")
    #with open("output.txt", "w") as file:
    #    for i in output_sequence_list:
    #        file.write(start + str(i) + end + '\n')


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
обратные им последовательности.

КАК ПОЛЬЗОВАТЬСЯ
Создайте в директории с программой директорию 'input', куда положите
Ваши файлы с расширением *.ab1
Введите начальную и конечную последовательность
Результат будет в файле output.fasta, список необработанных файлов
в error.txt
""")

subdir = 'input'
start = input('Введите начало последовательности: ').upper()
end = input('Введите конец последовательности: ').upper()
output_sequence_list = []
output_error_list = []
filename_list = listdir(subdir)

for filename in filename_list:
    path_file = path.join(subdir, filename)
    # Получаем 4 вида одной последовательности
    sequence = get_sequence(path_file)
    sequence_reverse = sequence[::-1]
    comp_sequence = make_complimentary(sequence)
    comp_sequence_reverse = comp_sequence[::-1]

    sequence_find = find_sequence(sequence, sequence_reverse,
                                  comp_sequence, comp_sequence_reverse,
                                  start, end)

    if sequence_find:
        # Записываем последовательность в список.
        output_sequence_list.append(SeqRecord(sequence_find))
        print(sequence_find)
    else:
        # Записываем название файла, в котором возникла ошибка.
        output_error_list.append(filename)

write_in_output_file(output_sequence_list)
write_in_error_file(output_error_list)
print('Готово!', len(output_sequence_list))

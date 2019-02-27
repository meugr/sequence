# Sequence

Программа для обработки последовательностей ДНК. Принимает на вход файлы последовательностей .ab1, на выход отдаёт файл в формате FASTA.

## Порядок работы
- принимает файлы *.ab1
- забирает из них последовательности ДНК
- ищет и вырезает последовательность между двумя фланкирующими последовательностями, включая их.
- возвращает файл fasta с искомыми последовательностями, где id=название исходного файла, descruption=некоторая информация о последовательности. Учитывает неоднозначные нуклеотиды.
- выводит в консоль полученные последовательности, по завершению работы сообщает общее количество входных файлов и количество удачно обработанных.

### Содержание description
1. Длина последовательности, включая фланкирующие.
2. Номера начального и конечного нуклеотидов искомой последовательности в родительской.
3. Тип последовательности (_sequence_  - прямая, _sequence_reverse_  - обратная, _comp_sequence_  - прямая комплементарная, _comp_sequence_reverse_  - обратная комплементарная).
4. Количество неоднозначных нуклеотидов в найденной последовательности.

## Установка
Для работы необходим `python 3.6` или новее, установленная библиотека `biopython` .

### Windows
- Скачайте и установите интерпретатор `python` с официального сайта [python.org](python.org)
- откройте `cmd`  и наберите там `python -m pip install biopython` 
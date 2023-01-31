import sys
import os
import subprocess
import cut_exon
import apply
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning) # FutureWarning 제거

FILE_MAX_SIZE = 2621440000
FILE_EXTENSION = ['.fastq', '.bam']

def main(argv):
    FILE_NAME = argv[1]
    
    path = os.path.join(os.getcwd(), FILE_NAME)
    file_size = os.path.getsize(path)
    root, extension = os.path.splitext(path)
    flag = 1

    if extension in FILE_EXTENSION:
        # print('정상적인 확장자입니다.')
        print('good extension.: {} file'.format(extension))
    else:
        # print('올바르지 않은 확장자입니다. 파일 입력을 다시 하세요.')
        print('incorrect extension. file must be .fastq or .bam')
        flag = 0
    
    if file_size <= FILE_MAX_SIZE:
        # print('파일 크기가 정상입니다.')
        print('file size is normal.')
    else :
        # print('파일 크기가 기준을 초과했습니다. 분석을 진행할 수 없습니다.')
        print('file size is over the range.')
        flag = 0
    
    size = file_size / (1024.0 * 1024.0 * 1000.0) ## GB로 출력
    print(f'{size:.2f} GB')
    
    if flag == 1:
        # make_consensus(FILE_NAME, extension)
        
        file_name = get_file_name(FILE_NAME, extension)
        f = open(os.path.join(os.getcwd(), 'exon_{}.txt'.format(file_name)), 'r')
        data = f.readline()
        apply.apply(data)


def make_consensus(file, ext):
    file_name = get_file_name(file, ext)
    if ext == '.fastq':
        subprocess.call(['./start_with_fastq.sh', file_name])
        cut_exon.cut(file_name)
    elif ext == '.bam':
        subprocess.call(['./start_with_bam.sh', file_name])
        cut_exon.cut(file_name)
    else:
        pass
    
    
def get_file_name(file, ext):
    if ext == '.fastq':
        file_name = file.rstrip('.fastq')
    elif ext == '.bam':
        file_name = file.rstrip('.bam')
    else:
        pass
    return file_name
    

if __name__ == "__main__":
    main(sys.argv)

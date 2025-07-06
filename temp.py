
    
import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np



#  файл и установка рабочей директории 
os.chdir(r"C:\Users\Admin\snp_project\PheGenI.csv")  
os.chdir(r'C:/Users/Admin/.spyder-py3')

# Загружаем CSV смотрим первые 5 строк и названия колонок
df = pd.read_csv("PheGenI.csv")
print(df.head())        
print(df.columns)         
    


# Преобразуем "Chromosome" в строку
df['Chromosome'] = df['Chromosome'].astype(str)

# Считаем количество SNP по каждой хромосоме
#chrom_counts = df['Chromosome'].value_counts().sort_index()

# Рисуем гистограмму
#plt.figure(figsize=(12, 6))
#chrom_counts.plot(kind='bar', color='cornflowerblue')
#plt.title("Количество SNP по хромосомам")
'''plt.xlabel("Хромосома")
plt.ylabel("Количество SNP")
plt.grid(axis='y')
plt.tight_layout()
plt.show()
  '''  

# Удалим странные значения, которые не являются номерами 1–22, X или Y
'''valid_chromosomes = [str(i) for i in range(1, 23)] + ['X', 'Y']
df = df[df['Chromosome'].isin(valid_chromosomes)]

# Перестроим график:
chrom_counts = df['Chromosome'].value_counts().sort_index()
chrom_counts.plot(kind='bar', figsize=(12,6), color='lightseagreen')
plt.title("Очищенное распределение SNP по хромосомам")
plt.xlabel("Хромосома")
plt.ylabel("Количество SNP")
plt.tight_layout()
plt.show()'''


valid_chromosomes = [str(i) for i in range(1, 23)] + ['X', 'Y']
df = df[df['Chromosome'].isin(valid_chromosomes)]

# Построим график с X и Y
chrom_counts = df['Chromosome'].value_counts().sort_index(key=lambda x: [int(i) if i.isdigit() else (23 if i == 'X' else 24 if i == 'Y' else 25) for i in x])
chrom_counts.plot(kind='bar', figsize=(12,6), color='cornflowerblue')
plt.title("SNPs by chromosome (including X and Y)")
plt.xlabel("Сhromosome")
plt.ylabel("Number of SNPs")
plt.tight_layout()
plt.show()

#TRAIT признаки
# Топ-10 признаков, с которыми связано больше всего SNP
trait_counts = df['Trait'].value_counts().head(10)



# Горизонтальный график
plt.figure(figsize=(10, 6))
trait_counts.plot(kind='barh', color='skyblue')
plt.title("Top 10 traits by SNP count")
plt.xlabel("Number of SNPs")
plt.ylabel("(Trait)")
plt.tight_layout()
plt.show()




print(df['P-Value'].head(10))

df['P-Value'] = pd.to_numeric(df['P-Value'], errors='coerce')  # пропускаем нечисловые
print(df['P-Value'].isna().sum(), "Missing values in P-Value column")
#statistic
print(df['P-Value'].describe())




# Заменим P-Value = 0 на очень маленькое число, чтобы не получить -inf
df['P-Value'] = df['P-Value'].replace(0, 1e-300)

# Преобразуем в -log10
neg_log_p = -np.log10(df['P-Value'])

#  график
plt.figure(figsize=(10, 6))
plt.hist(neg_log_p, bins=50, color='darkred')
plt.title("Distribution of SNP significance (-log10(P-Value))")
plt.xlabel("-log10(P-Value)")
plt.ylabel("Number of SNPs")
plt.tight_layout()
plt.show()


#какие именно SNP имеют P < 1e-10 (топ-10 самых значимых SNP)
top_snps = df[df['P-Value'] < 1e-10]
print(top_snps[['SNP rs', 'Gene', 'Trait', 'P-Value']].head(10))






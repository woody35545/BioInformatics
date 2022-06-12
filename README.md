# Predict non-small cell lung cancer using DNA sequence  
  

## [보유 Sequence Data]
#### index 1 - 1342 까지 비소세포폐암 환자 (1341 개)  
#### index 1342 - 2542 까지 일반인 (1201개)  
> 약 250명분의 Real Data와 원본 데이터에 노이즈를 추가한 Data Augmentation을 사용하여 총 2542개의 데이터로 구성하였음.

## [학습에 사용된 데이터]
> 환자 1200명, 일반인 1200명

2400개의 DNA Sequence Data를 한번에 학습하기에는 GPU Resource 관련 Issue가 발생하여서,
단계당 800명 데이터를 이용하여 3단계로 나눠서 학습하였음.(진행중)

### 1st train: 환자 400명 + 정상인 400명  
	index: (1 ~ 400) + (1343 ~ 1742)   
  
### 2nd train: 환자 400명 + 정상인 400명  
	index: (401 ~ 800) + (1743 ~ 2142)  
  
### 3rd train: 환자 400명 + 정상인 400명  
 	index: (801 ~ 1200) + (2143 ~ 2542)  

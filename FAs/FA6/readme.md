# Formative Assessment 6: Discrimination and Classification (LDA and QDA)

**Course:** APM1220 – Applied Multivariate Data Analysis

**## Joseph Raphael Espiritu ##**
**## Spike Patayon ##**

## 1. Project Description

This project investigates whether three varieties of wheat—Kama, Rosa, and Canadian—can be distinguished based on their geometrical characteristics using statistical discriminant analysis.

The analysis applies Linear Discriminant Analysis (LDA) and Quadratic Discriminant Analysis (QDA) to classify wheat varieties using seven continuous numerical predictors.

The main objective is to compare the predictive performance of both classification methods and determine which model provides a more suitable approach for identifying wheat varieties.

## 2. Dataset Source

**Dataset:** Seeds Dataset

**Source:** UCI Machine Learning Repository

**Link:** https://archive.ics.uci.edu/dataset/236/seeds

The dataset contains 210 observations, consisting of 70 observations for each wheat variety.

The seven numerical predictors are:

1. Area
2. Perimeter
3. Compactness
4. Kernel Length
5. Kernel Width
6. Asymmetry Coefficient
7. Kernel Groove Length

The response variable is Wheat Variety, consisting of three categorical classes: Kama, Rosa, and Canadian.

## 3. Methodology

The assessment is organized into five major parts:

- **Part 1 – Data Exploration and Preparation:** Dataset inspection, descriptive statistics, class distribution, missing values, and potential outliers.
- **Part 2 – Exploratory Analysis and Assumption Checking:** Graphical analysis, multivariate normality testing using the Henze–Zirkler test, and equality of covariance matrices using Box's M-test.
- **Part 3 – Model Development:** Stratified 70:30 training-testing partition, LDA and QDA model fitting, predicted class memberships, and posterior probabilities.
- **Part 4 – Model Evaluation and Comparison:** Confusion matrices, accuracy, misclassification rate, class-wise sensitivity and specificity, and stratified five-fold cross-validation.
- **Part 5 – Statistical Interpretation and Conclusion:** Interpretation of discriminating characteristics, classification boundaries, covariance structures, predictive performance, model recommendation, and limitations.

## 4. Technologies and Libraries

The analysis is implemented in Python using Google Colab and the following libraries:

- Pandas – Data manipulation and descriptive statistics
- NumPy – Numerical computations
- Matplotlib – Data visualization
- Seaborn – Statistical graphics
- SciPy – Statistical computations
- Pingouin – Multivariate normality and Box's M-test
- Scikit-learn – LDA, QDA, data partitioning, cross-validation, and model evaluation

## 5. Repository Contents

| File | Description |
|---|---|
| `FA6_Espiritu_Patayon.ipynb` | Python notebook containing the source code, statistical analysis, visualizations, and interpretations. |
| `FA6_Espiritu_Patayon.html` | Rendered HTML version of the completed statistical report. |
| `seeds_dataset.txt` | Original dataset obtained from the UCI Machine Learning Repository. |
| `README.md` | Project documentation and instructions for reproducing the analysis. |

## 6. Instructions for Reproducing the Results

1. Open [Google Colab](https://colab.research.google.com/).
2. Upload or open the provided `.ipynb` notebook.
3. Ensure that the required Python libraries are installed using the package installation cell.
4. Upload `seeds_dataset.txt` to the Colab runtime if prompted.
5. Execute the notebook cells sequentially, from the beginning to the end.
6. Review the generated statistical outputs, tables, visualizations, and interpretations.

A fixed random seed of 42 is used for reproducibility of the training-testing partition and cross-validation procedure.

The analysis uses the prescribed UCI Seeds Dataset without generating or simulating additional observations.

---

**Assessment:** Formative Assessment 6 – Discrimination and Classification (LDA and QDA)  
**Dataset Provider:** UCI Machine Learning Repository

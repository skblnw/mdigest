import numpy as np
from sklearn import svm
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score, roc_auc_score, roc_curve
import matplotlib.pyplot as plt

# Load your dataset
def load_data(filename):
    data = []
    with open(filename, 'r') as file:
        for line in file:
            parts = line.strip().split()
            label = int(parts[0])
            features = [float(x) for x in parts[1:]]
            data.append((features, label))
    return data

# Split dataset into features and labels
def split_features_labels(data):
    features = [x[0] for x in data]
    labels = [x[1] for x in data]
    return np.array(features), np.array(labels)

# Main function for SVM with Stratified 5-Fold cross-validation
def svm_stratified_fold_cross_validation(data_filename):
    # Load and prepare data
    data = load_data(data_filename)
    X, y = split_features_labels(data)

    # Define Stratified 5-Fold cross validation
    skf = StratifiedKFold(n_splits=5)

    # Initialize accuracy and AUC lists
    accuracies = []
    auc_scores = []

    # Perform Stratified 5-Fold cross-validation
    for train_index, test_index in skf.split(X, y):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]

        # Create SVM classifier with probability estimation
        clf = svm.SVC(kernel='linear', probability=True)

        # Train the model
        clf.fit(X_train, y_train)

        # Make predictions and compute accuracy
        y_pred = clf.predict(X_test)
        accuracy = accuracy_score(y_test, y_pred)
        accuracies.append(accuracy)

        # Compute AUC-ROC
        y_scores = clf.predict_proba(X_test)[:, 1]
        auc_score = roc_auc_score(y_test, y_scores)
        auc_scores.append(auc_score)

        print(f"Fold accuracy: {accuracy}, AUC: {auc_score}")

        # Plot ROC curve
        fpr, tpr, thresholds = roc_curve(y_test, y_scores)
        plt.plot(fpr, tpr, color='blue', label=f'ROC Curve (AUC = {auc:.2f})')
        plt.plot([0, 1], [0, 1], color='grey', linestyle='--')
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('Receiver Operating Characteristic (ROC) Curve')
        plt.legend()
        plt.show()

    # Print average accuracy and AUC
    print(f"Average Accuracy: {np.mean(accuracies)}")
    print(f"Average AUC: {np.mean(auc_scores)}")

# Example usage
svm_stratified_fold_cross_validation('gcc_submatrix.txt')

import numpy as np
from sklearn import svm
from sklearn.model_selection import KFold
from sklearn.metrics import accuracy_score

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

# Main function for SVM with 5-Fold cross-validation
def svm_fold_cross_validation(data_filename):
    # Load and prepare data
    data = load_data(data_filename)
    X, y = split_features_labels(data)

    # Define 5-Fold cross validation
    kf = KFold(n_splits=5)

    # Initialize accuracy list
    accuracies = []

    # Perform 5-Fold cross-validation
    for train_index, test_index in kf.split(X):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]

        # Create SVM classifier
        clf = svm.SVC(kernel='linear')

        # Train the model
        clf.fit(X_train, y_train)

        # Make predictions
        y_pred = clf.predict(X_test)

        # Calculate accuracy
        accuracy = accuracy_score(y_test, y_pred)
        accuracies.append(accuracy)
        print(f"Fold accuracy: {accuracy}")

    # Print average accuracy
    print(f"Average Accuracy: {np.mean(accuracies)}")

# Example usage
# svm_fold_cross_validation('your_dataset.txt')

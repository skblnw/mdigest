import numpy as np
from sklearn import svm
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score, roc_auc_score, roc_curve
import matplotlib.pyplot as plt
import joblib  # for saving the model

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
def svm_stratified_fold_cross_validation(data_filename, model_filename):
    # Load and prepare data
    data = load_data(data_filename)
    X, y = split_features_labels(data)

    # Define Stratified 5-Fold cross validation
    skf = StratifiedKFold(n_splits=5)

    # Initialize accuracy and AUC lists
    best_model = None
    best_accuracy = 0
    accuracies = []
    auc_scores = []

    # Perform Stratified 5-Fold cross-validation
    for train_index, test_index in skf.split(X, y):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]

        # Create SVM classifier with probability estimation
        model = svm.SVC(kernel='linear', probability=True)

        # Train the model
        model.fit(X_train, y_train)

        # Make predictions and compute accuracy
        y_pred = model.predict(X_test)
        accuracy = accuracy_score(y_test, y_pred)
        accuracies.append(accuracy)

        # Compute AUC-ROC
        y_scores = model.predict_proba(X_test)[:, 1]
        auc_score = roc_auc_score(y_test, y_scores)
        auc_scores.append(auc_score)

        print(f"Fold accuracy: {accuracy}, AUC: {auc_score}")

        # Plot ROC curve
        fpr, tpr, thresholds = roc_curve(y_test, y_scores)
        plt.plot(fpr, tpr, color='blue', label=f'ROC Curve (AUC = {auc_score:.2f})')
        plt.plot([0, 1], [0, 1], color='grey', linestyle='--')
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('Receiver Operating Characteristic (ROC) Curve')
        plt.legend()
        # plt.show()

        if accuracy > best_accuracy:
            best_accuracy = accuracy
            best_model = model

    # Print average accuracy and AUC
    print(f"Average Accuracy: {np.mean(accuracies)}")
    print(f"Average AUC: {np.mean(auc_scores)}")
    joblib.dump(best_model, model_filename)
    print(f"Best model saved with accuracy: {best_accuracy}")

# Function to test the saved model on a new dataset
def test_saved_model(model_filename, test_data_filename):
    test_data = load_data(test_data_filename)
    X_test, y_test = split_features_labels(test_data)
    print(y_test)

    loaded_model = joblib.load(model_filename)
    y_pred = loaded_model.predict(X_test)
    print(y_pred)
    test_accuracy = accuracy_score(y_test, y_pred)
    print(f"Test accuracy on new dataset: {test_accuracy}")

# Example usage
model_file = 'best_svm_model.pkl'
train_data_file = '/home/kevin/data/maap/md.a2/ILAKFLHWL-ila1.5men/hlaonly/analysis/mdigest/t1/gcc_submatrix.txt'
# train_data_file = 'gcc_submatrix.txt'
svm_stratified_fold_cross_validation(train_data_file, model_file)

test_data_file = '/home/kevin/data/maap/md.a2/EAAGIGILTV-mel8.7q9b/hlaonly/analysis/mdigest/t1/gcc_submatrix.txt'
# test_data_file = 't1/gcc_submatrix.txt'
test_saved_model(model_file, test_data_file)

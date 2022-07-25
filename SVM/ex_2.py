from sklearn.svm import SVC
classifier = SVC(kernel = 'linear')
training_points = [[1, 2], [1, 5], [2, 2], [7, 5], [9, 4], [8, 2]]
labels = [1, 1, 1, 0, 0, 0]
classifier.fit(training_points, labels) 
print(classifier.predict([[3, 2]]))

# classifier = SVC(C = 0.01)
# classifier = SVC(kernel = 'linear')
# classifier = SVC(kernel='rbf')
# classifier = SVC(kernel = "rbf", C = 2, gamma = 0.5)
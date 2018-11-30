def generate_data(n_samples, n_reps, desired_mean, desired_std, indices):
    np.random.seed(55)
    random_data = np.random.normal(loc=0.0, scale=desired_std, size=(n_reps, n_samples))
    new_data = pd.DataFrame(((random_data - np.mean(random_data)) * (desired_std/(np.std(random_data-np.mean(random_data)))) + desired_mean).transpose(), index=indices)
    return new_data

def take_isolated_data(data, corr_data, value_count):
    none = data.loc[value_count == 0]
    one = data.loc[value_count == 1]
    two = data.loc[value_count >= 2]
    
    #In order of 0, then 1, then 2 or more:
    tup_data = (none, one, two)
    tup_means = (np.min(data), np.mean(tup_data), np.mean(corr_data.loc[tup_data.index]))
    tup_stds = (np.std(corr_data.loc[tup_data.index], np.std(corr_data.loc[tup_data.index], np.std(two))))
    
    data.loc[tup_data.index] = generate_data(np.shape(tup_data)[0], np.shape(tup_data)[1], tup_means, tup_stds, tup_data.index)
    
    return data
    
    
    #Potential use of tuple unpacking of arguments.
    # new_data.loc[tup] = 


desired_stds = np.array([1, 2, 3, 1, 2])
desired_means = np.array([100, 300, 500, 700, 900])
names = pd.Series(['a', 'b', 'c', 'd', 'e'])

fake_data = generate_data(5, 3, desired_means, desired_stds, names)

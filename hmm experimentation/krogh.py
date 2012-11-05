import random

# A profile HMM as an accompaniment to Krogh et al's 1994 paper on them

class ProfileHiddenMarkovModel(object):
    def __init__(self, insertion_dist):
        self.insertion_dist = insertion_dist
        self.positions = list()

    def add_position(self, match_dist, mi, ii, di):
        self.positions.append(Position(self, match_dist, mi, ii, di))

    def link_positions(self, mm, md, id, im, dm, dd):
        '''Link up the last two positions'''
        self.positions[-2].add_link(self.positions[-1],
                                    mm, md, id, im, dm, dd)

    def add_final_position(self):
        self.final = FinalState(message='end state')
        self.positions[-1].complete(self.final)

    def walk(self):
        current_state = self.positions[0].match
        output = ''
        while current_state is not None:
            current_state.print_message()
            symbol = current_state.emit()
            if symbol is not None:
                output += symbol
            current_state = current_state.transition()
        print(output)

class Position(object):
    '''A position in the profile.'''
    def __init__(self, model, match_dist, mi, ii, di):
        '''mi: match state to insert state transition probability
        ii: insert state to insert state transition probability
        di: delete state to insert state transition probability'''
        self.match = State('match state', match_dist)
        self.insert = State('insertion state', model.insertion_dist)
        self.delete = State('delete state', None)

        # Link the states
        self.match.add_link(mi, self.insert)
        self.insert.add_link(ii, self.insert)
        self.delete.add_link(di, self.insert)

    def add_link(self, next_position, mm, md, id, im, dm, dd):
        '''mm: prob of transition from this position's match state
               to next position's match state

           id: prob of transition from this position's insert state to
               next position's delete state

           im: prob of transition from this position's insert state to next
               position's match state

           dm: prob of transition from this position's delete state to next
               position's match state
            
           dd: you can guess'''
        
        self.match.add_link(mm, next_position.match)
        self.match.add_link(md, next_position.delete)
        self.insert.add_link(id, next_position.delete)
        self.insert.add_link(im, next_position.match)
        self.delete.add_link(dm, next_position.match)
        self.delete.add_link(dd, next_position.delete)

    def complete(self, final_state):
        for state in (self.match, self.insert, self.delete):
            remainder = state.links[-1][0]
            state.add_link(remainder, final_state)
           

class State(object):
    def __init__(self, message, dist):
        '''Create a state with a given emission distribution. (Give None
        for a state that does not emit.) Emission distribution must be a
        list of (cutoff, symbol). During emission, a random number from
        0 to 1 will be generated, and the symbol with the lowest cutoff
        above the random number will be emitted.'''
        self.links = list()
        self.message = message
        self.dist = dist
        
    def add_link(self, prob, state):
        '''Add a transition to state "state" with transition probability
        "prob"'''
        # Mathematically, the transition probabilities are stored as a
        # cumulative distribution function. In Python terms, it's a list
        # of pairs (cutoff, state), and the difference between "cutoff"
        # and the previous cutoff is the probability of returning "state"
        if len(self.links) == 0:
            self.links.append((prob, state))
        else:
            prev_total = self.links[-1][0]
            cutoff = prev_total + prob
            if cutoff > 1.000001:
                raise ValueError('Transition probabilities sum to higher '
                                 + 'than one')
            self.links.append((cutoff, state))

    def print_message(self):
        # Just for fun while I'm watching the chain
        print(self.message)

    def emit(self):
        # No distribution stored means state does not emit symbols
        if self.dist is None:
            return None
        
        # Generate a symbol from the cumulative distribution function
        random_number = random.random()
        for cutoff, symbol in self.dist:
            if random_number < cutoff:
                return symbol
    
    def transition(self):
        # Retrieve a state from the cumulative distribution function
        random_number = random.random()
        for cutoff, state in self.links:
            if random_number < cutoff:
                return state

class FinalState(State):
    def __init__(self, message):
        State.__init__(self, message, None)

    def transition(self):
        return None

# Generate a hidden Markov model.

# Make a hidden markov model where insertions have a uniform distribution
# of emissions:
residues = ['A', 'L', 'I', 'V', 'P',
            'C', 'S', 'T', 'N', 'Q',
            'D', 'E', 'K', 'R', 'H',
            'Y', 'F', 'W', 'G', 'M']
example = ProfileHiddenMarkovModel([(.05*p, r) \
                                    for p, r in zip(range(1,21), residues)])

# Make the 0 position from Figure 1, which has no delete state and no
# emissions from the match state
example.add_position(None, .1, .5, 0)

# Make the first real position and link it up
example.add_position([(.5, 'A'), (1, 'V')], .1, .5, .1)
example.link_positions(mm=.9, md=0, id=.2, im=.3, dm=0, dd=0)

# Make the second real position and link it up
example.add_position([(.5, 'L'), (1, 'I')], mi=.1, ii=.5, di=.1)
example.link_positions(mm=.8, md=.1, id=.1, im=.4, dm=.8, dd=.1)

# And that's all I have the patience for, so tie it up
example.add_final_position()
